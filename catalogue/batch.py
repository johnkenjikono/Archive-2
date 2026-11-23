"""
Message Batches transport for steps.metadata_curation.CuratedLabeler: the same tiers and the same
requests as the synchronous path, sent through the Batch API at half price.

Safety rails, because this spends real money:
  - a hard budget: before a batch is submitted its cost is estimated (conservatively, ignoring the
    prompt cache) and nothing is sent if spent + estimate would pass the budget
  - a ledger (spend.jsonl) of what every finished batch actually cost, from each response's usage
  - submitted batch ids are saved (batches.json) before waiting on them, keyed on the exact requests,
    so an interrupted run picks up the same batches instead of paying for them twice

Batch requests cannot use server-side refusal fallbacks, so chunks that come back refused (Sonnet and
Opus), cut off at max_tokens, errored or expired are retried on the synchronous path, which has them.
A refused Haiku chunk is not retried: its records escalate to the Sonnet tier as usual.
"""
import datetime
import hashlib
import json
import os
import time
from pathlib import Path

try:
    from steps.metadata_curation import BATCH_DISCOUNT, PRICES, usage_cost
except ImportError:  # pragma: no cover
    from metadata_curation import BATCH_DISCOUNT, PRICES, usage_cost

MAX_REQUESTS_PER_BATCH = 10_000
MAX_BYTES_PER_BATCH = 200 * 1024 * 1024
# Output tokens per record used for budget checks only: deliberately high, thinking included.
# Measured 2026-11-21 per record: Haiku 157, Sonnet 74 per run, Opus 139; Fable 73 (low effort) and 151 (high)
# for labels, 229-264 per judged record.
OUTPUT_PER_RECORD = {"claude-haiku-5-5": 300, "claude-sonnet-5-5": 450, "claude-opus-5-5": 900,
                     "claude-fable-5-2": 400}
CHARS_PER_TOKEN = 3.0  # conservative for JSON-heavy prompts


class BudgetExceeded(RuntimeError):
    pass


def _transient_errors():
    """API errors worth waiting out while polling (dropped connections, timeouts, overload, throttling)."""
    import anthropic
    return (anthropic.APIConnectionError, anthropic.APITimeoutError, anthropic.InternalServerError,
            anthropic.RateLimitError)


def _now():
    return datetime.datetime.now(datetime.timezone.utc).isoformat(timespec="seconds")


class Ledger:
    """Append-only record of money spent: one JSON line per finished batch (per model) or sync retry."""

    def __init__(self, path):
        self.path = Path(path)
        self.entries = []
        if self.path.exists():
            with open(self.path, encoding="utf-8") as f:
                self.entries = [json.loads(line) for line in f if line.strip()]

    @property
    def total(self):
        return sum(e["usd"] for e in self.entries)

    def add(self, **entry):
        entry = {"time": _now(), **entry}
        self.entries.append(entry)
        self.path.parent.mkdir(parents=True, exist_ok=True)
        with open(self.path, "a", encoding="utf-8") as f:
            f.write(json.dumps(entry, sort_keys=True) + "\n")
        return entry


def estimate_usd(requests):
    """Upper-end USD for a list of batch requests (input from prompt size, output per record)."""
    usd = 0.0
    for request in requests:
        params = request["params"]
        model = params["model"]
        price_in, price_out = PRICES[model]
        tokens_in = len(json.dumps(params)) / CHARS_PER_TOKEN
        records = params["messages"][0]["content"].count("\nRecord ") + 1
        usd += (tokens_in * price_in + records * OUTPUT_PER_RECORD.get(model, 900) * price_out) / 1e6
    return usd * BATCH_DISCOUNT


class BatchTransport:
    def __init__(self, workdir, budget_usd, client=None, poll_seconds=60, sleep=time.sleep, log=print,
                 stage_prefix="", state_name="batches.json"):
        self.workdir = Path(workdir)
        self.budget = float(budget_usd)
        self.client = client
        self.poll_seconds = poll_seconds
        self.sleep = sleep
        self.log = log
        self.stage_prefix = stage_prefix
        self.ledger = Ledger(self.workdir / "spend.jsonl")  # shared by every command in the work directory
        self.state_path = self.workdir / state_name  # one per command, so concurrent commands never collide
        self.state = self._read_state()

    # --- state -----------------------------------------------------------------------------------
    def _read_state(self):
        try:
            with open(self.state_path, encoding="utf-8") as f:
                return json.load(f)
        except (FileNotFoundError, json.JSONDecodeError):
            return {}

    def _save_state(self):
        self.state = {**self._read_state(), **self.state}  # keep stages another process saved meanwhile
        self.state_path.parent.mkdir(parents=True, exist_ok=True)
        tmp = Path(f"{self.state_path}.tmp{os.getpid()}")
        with open(tmp, "w", encoding="utf-8") as f:
            json.dump(self.state, f, indent=1, sort_keys=True)
        os.replace(tmp, self.state_path)

    def _client(self):
        if self.client is None:
            import anthropic  # only needed when a model tier runs
            self.client = anthropic.Anthropic()
        return self.client

    def _patiently(self, what, call, attempts=40):
        """Read-only calls (status, results) retried through network trouble for up to about 3 hours.
        Batch creation is deliberately not retried here: a lost reply could otherwise submit it twice."""
        errors = _transient_errors()
        for attempt in range(attempts):
            try:
                return call()
            except errors as e:
                if attempt == attempts - 1:
                    raise
                delay = min(30 * (attempt + 1), 300)
                self.log(f"  {what}: {type(e).__name__}; retrying in {delay} s")
                self.sleep(delay)

    # --- transport interface ---------------------------------------------------------------------
    def call_many(self, labeler, jobs):
        """jobs: [(tier, {rid: fields}, disputes)]. Returns one {rid: answer} dict per job."""
        chunks = []
        for j, (tier, records, disputes) in enumerate(jobs):
            size = labeler.tiers[tier]["chunk"]
            ids = list(records)
            for start in range(0, len(ids), size):
                part = {i: records[i] for i in ids[start:start + size]}
                chunks.append({"id": f"j{j}c{start // size}", "job": j, "tier": tier, "records": part,
                               "disputes": disputes})
        results = [{} for _ in jobs]
        if not chunks:
            return results
        name = self.stage_prefix + "+".join(sorted({c["tier"] for c in chunks}))
        failed = self._send(labeler, chunks, name, results)
        if failed:
            self._retry(labeler, failed, results, name)
        return results

    def _send(self, labeler, chunks, name, results):
        """Run chunks as one batch stage; fill results and return the chunks that need a retry."""
        requests = [{"custom_id": c["id"],
                     "params": labeler.build_request(c["tier"], c["records"], c["disputes"], batch=True)}
                    for c in chunks]
        messages = self._run_stage(name, requests, labeler)
        failed = []
        for chunk in chunks:
            kind, message = messages.get(chunk["id"], ("missing", None))
            status = kind
            if kind == "succeeded":
                answers, status = labeler.parse_response(chunk["tier"], chunk["records"], message)
                if status == "ok":
                    results[chunk["job"]].update(answers)
                    continue
            chunk["status"] = status
            failed.append(chunk)
        return failed

    def _retry(self, labeler, failed, results, name):
        counts = {}
        for chunk in failed:
            counts[chunk["status"]] = counts.get(chunk["status"], 0) + 1
        self.log(f"  {name}: {len(failed)} chunk(s) to retry ({', '.join(f'{k} {v}' for k, v in sorted(counts.items()))})")
        retry = [c for c in failed if not (c["status"] == "refusal" and not labeler.tiers[c["tier"]].get("fallbacks"))]
        # Whole batches that errored or expired go back as one more batch; a few stragglers go synchronously.
        resend = [c for c in retry if c["status"] in ("errored", "expired", "canceled", "missing")]
        if len(resend) > 50:
            for c in resend:
                c["id"] = f"{c['id']}r"
            still = self._send(labeler, resend, f"{name} resend", results)
            retry = [c for c in retry if c not in resend] + still
        if not retry:
            return
        requests = [{"params": labeler.build_request(c["tier"], c["records"], c["disputes"])} for c in retry]
        estimate = estimate_usd(requests) / BATCH_DISCOUNT  # the synchronous path pays full price
        self._check_budget(estimate, f"{name} synchronous retries")
        before = labeler.cost()
        for chunk in retry:
            results[chunk["job"]].update(labeler._call_chunk(chunk["tier"], chunk["records"], chunk["disputes"]))
        spent = labeler.cost() - before
        self.ledger.add(stage=f"{name} retries", kind="sync", model="mixed", usd=round(spent, 9),
                        requests=len(retry))

    # --- batches ---------------------------------------------------------------------------------
    def _check_budget(self, estimate, what):
        self.ledger = Ledger(self.ledger.path)  # count what other commands have spent since this one started
        if self.ledger.total + estimate > self.budget:
            raise BudgetExceeded(f"{what} would cost up to ${estimate:.2f}; ${self.ledger.total:.2f} of the "
                                 f"${self.budget:.2f} budget is already spent. Raise --budget to continue.")

    @staticmethod
    def _groups(requests):
        group, size = [], 0
        for request in requests:
            n = len(json.dumps(request))
            if group and (len(group) >= MAX_REQUESTS_PER_BATCH or size + n > MAX_BYTES_PER_BATCH):
                yield group
                group, size = [], 0
            group.append(request)
            size += n
        if group:
            yield group

    def _run_stage(self, name, requests, labeler):
        """Submit (or resume) the batches for these exact requests; returns {custom_id: (type, message)}."""
        key = hashlib.sha256(json.dumps(requests, sort_keys=True).encode()).hexdigest()[:24]
        groups = list(self._groups(requests))
        stage = self.state.get(key)
        if stage is None:
            estimate = estimate_usd(requests)
            self._check_budget(estimate, f"{name} ({len(requests)} requests)")
            stage = self.state[key] = {"name": name, "batches": [], "groups": len(groups), "requests": len(requests),
                                       "estimate_usd": round(estimate, 4), "submitted": _now(), "accounted": []}
            self._save_state()
            self.log(f"  {name}: submitting {len(requests)} request(s) in {len(groups)} batch(es), "
                     f"estimated at most ${estimate:.2f}")
        else:
            self.log(f"  {name}: resuming {len(stage['batches'])} saved batch(es)")
        for group in groups[len(stage["batches"]):]:  # a run interrupted mid-submission sends only the rest
            batch = self._client().messages.batches.create(requests=group)
            stage["batches"].append(batch.id)
            self._save_state()
        out = {}
        for batch_id in stage["batches"]:
            self._wait(batch_id, name)
            spend = {}
            fetched = self._patiently(f"{name} results",
                                      lambda: list(self._client().messages.batches.results(batch_id)))
            for result in fetched:
                outcome = result.result
                message = getattr(outcome, "message", None)
                out[result.custom_id] = (outcome.type, message)
                if outcome.type == "succeeded" and message is not None:
                    model = getattr(message, "model", None) or "unknown"
                    labeler.record_usage(model, message.usage, batch=True)
                    s = spend.setdefault(model, {"usd": 0.0, "requests": 0, "input_tokens": 0, "output_tokens": 0,
                                                 "cache_read_input_tokens": 0, "cache_creation_input_tokens": 0})
                    s["usd"] += usage_cost(model, message.usage, batch=True)
                    s["requests"] += 1
                    for field in ("input_tokens", "output_tokens", "cache_read_input_tokens",
                                  "cache_creation_input_tokens"):
                        s[field] += getattr(message.usage, field, 0) or 0
            if batch_id not in stage["accounted"]:  # each batch is charged to the ledger once
                for model, s in spend.items():
                    self.ledger.add(stage=name, kind="batch", batch=batch_id, model=model, **{
                        k: (round(v, 9) if k == "usd" else v) for k, v in s.items()})
                stage["accounted"].append(batch_id)
                stage["finished"] = _now()
                self._save_state()
        return out

    def _wait(self, batch_id, name):
        last_log = 0.0
        while True:
            batch = self._patiently(f"{name} status", lambda: self._client().messages.batches.retrieve(batch_id))
            if batch.processing_status == "ended":
                counts = batch.request_counts
                self.log(f"  {name}: batch {batch_id} ended (succeeded {counts.succeeded}, errored {counts.errored}, "
                         f"expired {counts.expired}, canceled {counts.canceled})")
                return batch
            now = time.monotonic()
            if now - last_log > 600:
                counts = batch.request_counts
                self.log(f"  {name}: batch {batch_id} {batch.processing_status}, {counts.processing} processing, "
                         f"{counts.succeeded} succeeded")
                last_log = now
            self.sleep(self.poll_seconds)
