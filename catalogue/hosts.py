"""
Host names and taxids resolved through NCBI Taxonomy, then placed in the curation vocabulary.

A BioSample host field names the host in many ways ("Homo sapiens", "human", "Gallus gallus domesticus",
"Oryza sativa L.", "Ixodes scapularis", "pig; slaughterhouse"). For each distinct value this tries the
whole value, its first part before ; , or /, a name in parentheses, and the leading binomial, and keeps
the first candidate NCBI recognises as a eukaryote (hosts of bacteria are eukaryotes; this also stops
bacterial genus names that happen to match a word). The vocabulary host is then the most specific
vocabulary entry on that taxon's lineage: Ixodes scapularis -> tick (Ixodida), Sus scrofa domesticus ->
pig, Acyrthosiphon pisum -> insect, Oryza sativa Japonica Group -> rice.
"""
import re

try:
    from steps.metadata_curation import normalize
except ImportError:  # pragma: no cover - run from another directory
    from metadata_curation import normalize

EUKARYOTA = 2759


def candidates(value):
    """Name strings worth looking up for one host value, most literal first."""
    text = normalize(value)
    if not text:
        return []
    out = [text]
    first = normalize(re.split(r"[;,/]", text)[0])
    out.append(first)
    inner = re.search(r"\(([^)]*)\)", text)
    if inner:
        out.append(normalize(inner.group(1)))
    outside = normalize(re.sub(r"\([^)]*\)", " ", first))
    out.append(outside)
    words = outside.split()
    if len(words) > 2:
        out.append(" ".join(words[:2]))
    seen, ordered = set(), []
    for c in out:
        if c and c not in seen and len(c) > 2 and not c.isdigit():
            seen.add(c)
            ordered.append(c)
    return ordered


class HostResolver:
    """taxa: {query: record} from catalogue.ncbi.resolve_taxa over every candidate and host taxid."""

    def __init__(self, vocab, taxa):
        self.vocab = vocab
        self.taxa = taxa
        self.vocab_taxid = {int(h["taxid"]): key for key, h in vocab.hosts.items() if h.get("taxid")}

    @staticmethod
    def queries(host_values=(), host_taxids=()):
        """Every lookup the resolver may need, for resolve_taxa."""
        names = {c for v in host_values for c in candidates(v)}
        ids = {str(t).strip() for t in host_taxids if str(t).strip().isdigit()}
        return names | ids

    def _eukaryote(self, record, query=None):
        """The record if it is a single eukaryote taxon. A one-word query that matched a genus or higher
        by its scientific name is refused: English words such as "cancer" or "lens" are also genus names."""
        if not record or "ambiguous" in record:
            return None
        lineage = [int(t) for t in record.get("parents") or []] + [int(record["tax_id"])]
        if EUKARYOTA not in lineage:
            return None
        if (query and " " not in query and normalize(record.get("name", "")) == query
                and record.get("rank") not in ("species", "subspecies", "varietas", "forma", "strain")):
            return None
        return record

    def vocab_key(self, record):
        """Most specific vocabulary host on the lineage of a taxonomy record, or None."""
        if not record:
            return None
        lineage = [int(t) for t in record.get("parents") or []] + [int(record["tax_id"])]
        for taxid in reversed(lineage):
            if taxid in self.vocab_taxid:
                return self.vocab_taxid[taxid]
        return None

    def lookup_name(self, value):
        """The NCBI record for a host value, or None."""
        for candidate in candidates(value):
            record = self._eukaryote(self.taxa.get(candidate), candidate)
            if record:
                return record
        return None

    def lookup_taxid(self, taxid):
        return self._eukaryote(self.taxa.get(str(taxid).strip()))

    def vocab_key_for_name(self, value):
        return self.vocab_key(self.lookup_name(value))

    def vocab_key_for_taxid(self, taxid):
        return self.vocab_key(self.lookup_taxid(taxid))

    def exact(self, attributes, fields=("host",)):
        """(taxid, scientific name) of the host a BioSample names: host_taxid first, then the host fields."""
        record = None
        if str(attributes.get("host_taxid", "")).strip():
            record = self.lookup_taxid(attributes["host_taxid"])
        for field in fields:
            if record is None and attributes.get(field):
                record = self.lookup_name(attributes[field])
        return (record["tax_id"], record["name"]) if record else ("", "")
