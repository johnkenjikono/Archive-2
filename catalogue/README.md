# Isolate metadata catalogue

Curated host and isolation-source metadata for every bacterial genome NCBI holds at complete-genome or
chromosome level (108,400 GenBank assemblies on 2026-11-20, 8,200 of them metagenome-assembled), with every
value tied to the BioSample words it came from. It runs the pipeline's tiered curation
(`steps/metadata_curation.py`: rules, then Haiku, two agreeing Sonnet runs and Opus) over all of NCBI instead
of one species, so the labels are comparable across species and the pipeline's `--niche` step can reuse them.

## Running it

From the repository root (the work directory is created as `./catalogue_work`, gitignored, or wherever
`--workdir` points):

```bash
python -m catalogue.build harvest                       # NCBI assembly reports + BioSample attributes (free, ~2 min)
python -m catalogue.build prepare                       # taxonomy, distinct records, rules, cost estimate (free)
python -m catalogue.build run --budget 10 --sample 500  # pilot: 500 random open records through the Batch API
python -m catalogue.build run --budget 150              # every record the rules and the cache leave open
python -m catalogue.build verify --budget 30 --sample 500  # pilot: Claude Fable 5.2 cross-checks 500 records
python -m catalogue.build verify --budget 200 --effort high  # Fable cross-checks every model label, judges disputes
python -m catalogue.build rejudge --budget 100 --sample 300  # Opus 5.5 judges 300 of those disputes again
python -m catalogue.build assemble                      # the catalogue table, species summary and manifest
python -m catalogue.build status                        # what has run and what it has cost
```

It refuses to start if `steps/` would be loaded from a different checkout of the repository (possible
because `steps/` has no `__init__.py`, e.g. when a worktree's code is run from the main checkout's
folder). Only `run`, `verify` and `rejudge` call the Claude API.
It reads `ANTHROPIC_API_KEY` from the environment or from a `.env` file in the current directory or the
repository root; `.env` is gitignored. Everything is resumable: the harvest saves its page token, taxonomy
answers and model answers are cached, and submitted batch ids are saved before the run waits on them.

## What each stage does

| Stage | What happens | Cost |
|---|---|---|
| `harvest` | Pages through the NCBI Datasets v2 genome reports for taxon 2 (Bacteria), GenBank, complete genome + chromosome, keeping the accession, organism, assembly facts, CheckM, NCBI's taxonomy check and every BioSample attribute (`ncbi/pages/*.jsonl.gz`). Then lists the metagenome-assembled genomes among them | free |
| `prepare` | Looks up every host name, host taxid and organism taxid in NCBI Taxonomy (about 30,000 lookups the first time, cached in `taxa_cache.json`), groups genomes into distinct BioSample records, applies the rules tier and estimates what the models would cost (`reports/prepare_report.md`) | free |
| `run` | Sends the records the rules leave open through the model tiers using the Message Batches API at half price. A hard `--budget` (USD, for everything spent in the work directory) is checked before each batch, against a deliberately high estimate. Actual spend from each response's usage goes to `spend.jsonl`. Refused, errored or cut-off chunks are retried on the synchronous path, which has server-side refusal fallbacks | see below |
| `verify` | Claude Fable 5.2 labels the same records independently and a Fable judge settles disagreements (see "Cross-check" below). Same budget, ledger and resume rules as `run`; its batches are tracked in `crosscheck_batches.json` | see below |
| `assemble` | Labels every genome from the rules, the cached model answers and the cross-check's corrections, without calling any model (records no model has seen yet come out as `not curated`), and writes the release files | free |

### Which BioSample fields are read

The pipeline's fields (`host`, `host_disease`, `isolation_source`, `sample_type`, `source_type`,
`env_broad_scale`, `env_local_scale`, `env_medium`, `host_tissue_sampled`, `host_body_site`, `host_taxid`)
plus non-standard ones that, across the harvested genomes, carry host or source information:
`host_common_name`, `isolation_site`, `source`, `isolation_source_category`, `host_body_habitat`,
`plant_body_site`, `tissue`, `food_source`, `habitat`, `environment`, `metagenome_source` and
`coll_site_geo_feat`. Names are matched case-insensitively with spaces and hyphens read as `_`, and ENA's
names for the MIxS environment fields ("local environmental context", "environmental medium",
"broad-scale environmental context") count as `env_local_scale`, `env_medium` and `env_broad_scale`.
A genome's distinct record is the combination of these fields' values; genomes sharing a record share its
labels, so the models see each combination once.

### Hosts

Besides the vocabulary's own synonyms, the rules resolve any host name or taxid through NCBI Taxonomy and
place it on the most specific vocabulary host along its lineage: *Ixodes scapularis* becomes `tick`
(Ixodida), *Sus scrofa domesticus* `pig`, *Acyrthosiphon pisum* `insect`, *Oryza sativa* L. `rice`. Only
eukaryotic taxa count, and a one-word name that matched a genus by its scientific name is refused (English
words such as "cancer" are also genus names). The exact taxon named is kept in `host_ncbi_taxid` and
`host_ncbi_name`.

When the rules cannot read a record's sample description (it goes to the models anyway), they settle its host
only from a host taxid or a host field naming a particular host. A vague host ("animal", "plant", "soil") or a
host word in another field (`source_type` "animal") is left to the model reading the whole record, which names
the pig in "Rectum of pig". Fable's cross-check found the rules losing these cases to the record's own text.

A record that names a living host but has no source field at all gets the source "host sample, site not
stated" (or "plant, part not stated"), with the host field as its evidence. Its source category (human-,
animal- or plant-associated) is then known even though the body site is not. This applies to the catalogue
only; the pipeline's `--curated-metadata` leaves such sources unknown.

## Cross-check with Claude Fable 5.2 (`verify`)

The tiers above never check an answer Haiku is certain of, and that is most of them. `verify` has Claude
Fable 5.2, Anthropic's most capable model, label every record the curation sent to the models again,
independently: the same prompt, vocabulary, answer schema and quote check, without seeing the first answers.
Each part (host, source) is then compared with the curation's answer:

| Outcome | Meaning | What happens |
|---|---|---|
| `agrees` | The same host, or the same source term | Kept |
| `same niche` | (Source) a different term in the same niche | Kept as curated |
| `judge kept curation`, `judge took Fable's`, `judge chose another`, `judge: unknown` | The two disagree on a part the models decided | A judge (Fable at high effort) sees the record and both answers, in an order that alternates between records, through the curation's "two reviewers disagreed" prompt. Its answer replaces the curated one, or the value becomes unknown when it cannot back one with a quote |
| `rule disputed` | Fable disagrees with a rule's answer | Reported, not judged: a wrong rule is fixed in the rules or the vocabulary. What remains is mostly a host field that the source text contradicts (host *Bos taurus* on "beef carcass") |
| `no valid answer` | Fable's quote is not in the record | Kept as curated, unchecked |

`--effort` sets Fable's effort for its own labels (default `low`); the judge always runs at `high`. A run reports
the outcomes by curation tier, how often Fable's labels change between efforts, the cost per record and a
projection for every open record (`reports/verify_*.json`), and lists every disagreement with both answers
(`reports/verify_*_disagreements.tsv`). Labels and decisions are cached in `crosscheck_cache.json` per effort;
`assemble` uses the highest effort available for each record.

Full pass on 2026-11-22 (Fable at high effort, all 17,948 records the curation sent to the models, 44,511
genomes):

| Part | Fable agrees | Judge took Fable's | Judge kept curation | Judge: unknown | Rule disputed |
|---|---|---|---|---|---|
| Host | 17,236 (96%) | 422 | 141 | 86 | 63 |
| Source | 15,361 (86%, 769 of them same niche) | 2,125 (+14 another term) | 330 | 50 | 68 |

The judge changed the source of 9.5% of the records Haiku had settled, 16.5% of those settled by two Sonnet
runs and 17% of Opus's. It also found a source for 95 of the 126 records the curation had withheld. Across the
catalogue this moved 6,119 genomes (5.6%) to another niche, most often from unknown, laboratory culture or a
built environment (a hospital's name given as the source) to a host sample. A read of 30 random source changes
found about 21 clear fixes, 8 equally defensible and 1 worse. The judge took Fable's answer 80–81% of the
time whether that answer was shown first or second.

The pilots (`verify --sample 500`, $2.21 at low effort, $2.77 at high) agreed with these rates. High and low
effort differ on 8% of Fable's source labels, high more often right.

### A second judge (`rejudge`)

The judge is also Fable, so it may favour Fable's answers. `rejudge` tests this by having Opus 5.5 (`--model`)
judge a random sample of the settled disputes (`--sample`, default 300). Opus runs at the same effort, on the same
prompts, with the answers in the same order, so the model is the only difference between the two judges. The
sample leaves out records whose curated answer came from Opus, so Opus is a third party to both answers. Its
decisions are cached apart and never change the catalogue. `reports/rejudge_*.json` gives how often the two judges
agree: overall, per part, per curation tier, and by whether Fable was certain. `reports/rejudge_*.tsv` lists every
part with both decisions.

A Fable judge that favoured Fable would often take Fable's answer where Opus keeps the curation's, and rarely the
reverse. An exact McNemar test compares those two counts.

Run on 2026-11-23 for $1.03. Opus 5.5 at high effort judged 300 records with 331 disputed parts, sampled from the
2,421 of 2,887 judged records that have no curated answer from Opus:

| Fable's judge | Opus 5.5 | Parts |
|---|---|---|
| took Fable's answer | took Fable's answer | 246 |
| kept the curation's | kept the curation's | 30 |
| took Fable's answer | kept the curation's | 38 |
| kept the curation's | took Fable's answer | 13 |
| either judge chose a third answer | | 4 |

The two judges made the same decision on 276 parts (83%). Opus backed 246 of the 285 changes that Fable's judge
made toward Fable's answer (86%). Fable's judge does take Fable's answer more often than Opus does: 38 against 13
(McNemar p = 0.0006; host 12 against 3, source 26 against 10).

The difference comes entirely from answers Fable had marked as not certain:

| Fable's own answer | Parts | Same decision | Only Fable's judge took Fable's answer | Only Opus took it |
|---|---|---|---|---|
| Certain | 112 | 96% | 2 | 2 |
| Not certain | 219 | 77% | 36 | 11 |

By curation tier, the judges agree on 87% of Haiku-settled parts and 79% of parts settled by two Sonnet runs.

The 51 parts where only one judge took Fable's answer were then read by Claude Opus 5.5. That is the second
judge's own model, so this reading is not ground truth either:

- **Fable's answer better:** 16
- **The curation's better:** 6
- **Equally defensible:** 28
- **Unclear:** 1

The 38 parts where only Fable's judge took Fable's answer are the ones that would show self-preference. In those,
Fable's judge was right in 11 and Opus in 6.

Several of the 38 are records that give only ontology codes, such as env_medium `envo_00002262`. Fable reads the
code, while Opus answers unknown. Checked against EBI OLS, the codes support Fable's reading in 6 of 8 such
parts; `envo_00002262`, for example, is clay soil. So Fable's judge does lean toward Fable's answers, but its extra
changes are more often right than wrong, and the catalogue keeps its decisions.

## Cost

Measured on 2026-11-21 (vocabulary 2026.11.2): of 30,249 distinct records the rules settle 11,018 (63,883
genomes, of which 10,395 have no source or host field at all), and 17,954 records covering 41,362 genomes go to
the models. The full `run` took 12 minutes and cost $12.59 at Batch API prices (half of list):

| Tier | Records | Cost |
|---|---|---|
| Haiku on every open record | 17,954 (11,540 accepted) | $0.86 |
| Two Sonnet runs on the rest | 6,414 | $7.36 |
| Opus where the Sonnet runs disagree | 1,715 | $4.38 |

`prepare` estimates these costs beforehand from token counts (it expected $24.65 for this run), and a pilot
(`run --sample`) measures the real rates and projects the full run (`reports/run_pilot*.json`). The same holds
for `verify --sample`.

The cross-check (`verify --effort high`, 2026-11-22) took 6 minutes and cost $77.79:

| Step | Records | Cost |
|---|---|---|
| Fable labels at high effort | 17,451 (plus 497 from the pilot) | $63.30 |
| Fable judge at high effort | 2,821 with a model-decided disagreement | $14.49 |

Its budget pre-check assumes 400 output tokens per record. That is three times the 131 measured per label and
twice the 193 per judged record, so it needs `--budget 200` on a work directory that has also paid for `run`.

The second judge (`rejudge`, Opus 5.5 on 300 records) cost $1.03 against a pre-check estimate of $3.14. That
brings the work directory's spending to $96.91.

## Outputs (`catalogue_work/release/`)

`isolate_catalogue.tsv.gz`, one row per assembly:

| Columns | Meaning |
|---|---|
| `accession`, `paired_accession`, `biosample`, `bioproject` | GenBank (GCA) and RefSeq (GCF) accessions, BioSample, BioProject |
| `organism_name`, `tax_id`, `species`, `species_taxid`, `genus` … `phylum` | NCBI organism and its lineage |
| `strain`, `isolate`, `serovar` | As submitted |
| `assembly_level`, `release_date`, `is_mag`, `atypical`, `type_material`, `checkm_completeness`, `checkm_contamination`, `taxonomy_check` | Assembly facts and NCBI quality flags |
| `host`, `host_taxid` | Vocabulary host (a species such as `human` or a group such as `insect`) and its NCBI taxid |
| `host_ncbi_name`, `host_ncbi_taxid` | The exact taxon the record names, when NCBI recognises it |
| `source_category`, `mixs_extension`, `niche`, `source_term`, `source_term_id`, `body_site` / `food` / `environment` | Where the isolate came from: one of 18 source categories, one of 41 niches, and an ENVO, UBERON, PO or FoodOn term |
| `*_evidence_field`, `*_quote`, `*_tier`, `*_model`, `*_confidence` | For host and source: the field and the exact words each value rests on, the tier that decided it (`rule`, `haiku`, `sonnet x2`, `opus`, `fable judge`, `withheld`, `withheld by judge`, `not curated`, `no evidence`), the model, and the tier's confidence (`high`, `medium`, `low`, `none`) |
| `host_crosscheck`, `source_crosscheck`, `fable_host`, `fable_source_term`, `crosscheck_model` | The cross-check outcome per part (see "Cross-check"; `rules only` where the rules decided the whole record, `not checked` before `verify`), Fable's own answers, and the Fable model and effort |
| `country`, `country_insdc`, `collection_year` | Parsed from `geo_loc_name` and `collection_date`; the country checked against the INSDC list |
| `raw_*` | The BioSample values the curation read |
| `vocabulary_version` | Version of `steps/metadata_vocabulary.json` |

`species_summary.tsv`, one row per species: genomes, isolates, complete genomes, type material, genomes with a
known source, number of source categories and their Shannon entropy, top categories, niches and hosts,
countries, first and last collection year, and BioProjects. It doubles as a ranking of species by ecological
breadth for choosing pipeline species.

`manifest.json`: the NCBI query and harvest dates, counts, vocabulary version, models per tier, code commit and
total spend.

## Limits

- Only BioSample text is read. Genomes whose BioSample says nothing about the source (about one in ten here)
  stay `unknown`; the linked papers are not read.
- The confidence values are tiers of the method, not measured error rates. The cross-check measures agreement
  between models, not accuracy, and models can share a mistake. The Fable judge showed no order bias, but it
  takes Fable's answer more often than an Opus 5.5 judge does (38 against 13 parts in a sample of 331). The
  extra changes looked more often right than wrong, though that check was a model's reading too. An error rate
  needs a hand-labelled random sample.
- 151 records (332 genomes) give only ontology codes in their source fields, such as env_medium
  `envo_00002262`. The rules accept a bare code only when it is a vocabulary term's own id, so most of these
  records go to the models. The curation answered unknown for 119 of them, and Fable's judge has since given
  110 a term. Mapping each code through the ontology to its nearest vocabulary term (clay soil is a soil) would
  settle them by rule, for free.
- The rules take an environment field over a named host in 123 records (285 genomes), such as `env_medium`
  "seawater" on a marine sponge. Some of these are right (soil sampled around a named plant) and some are not.
- The vocabulary has no host group for protists or algae (ciliates, dinoflagellates, diatoms, amoebae). Their
  endosymbionts get an unknown host, or a wrong "animal" or "plant" one where the judge did not see them.
- For metagenome-assembled genomes the labels describe the metagenome sample, not an isolate.
- The vocabulary is fixed and coarse for some environments; free-text detail stays in the `raw_*` columns and
  the quotes.
