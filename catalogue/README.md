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
python -m catalogue.build prepare                       # taxonomy + ontology lookups, rules, cost estimate (free)
python -m catalogue.build run --budget 10 --sample 500  # pilot: 500 random open records through the Batch API
python -m catalogue.build run --budget 150              # every record the rules and the cache leave open
python -m catalogue.build verify --budget 30 --sample 500  # pilot: Claude Fable 5.2 cross-checks 500 records
python -m catalogue.build verify --budget 200 --effort high  # Fable cross-checks every model label, judges disputes
python -m catalogue.build rejudge --budget 100 --sample 300  # Opus 5.5 judges 300 of those disputes again
python -m catalogue.build audit-rules --budget 225  # Fable checks every record the rules settle; reported only
python -m catalogue.build papers-pilot --budget 146 --sample 200  # Fable reads 200 linked papers
python -m catalogue.build papers-fill --budget 600 --dry-run  # price reading unknown-field genomes' papers
python -m catalogue.build structure-pilot --budget 175 --species "Prochlorococcus marinus"  # published lineages
python -m catalogue.build assemble                      # the catalogue table, species summary and manifest
python -m catalogue.build status                        # what has run and what it has cost
```

It refuses to start if `steps/` would be loaded from a different checkout of the repository (possible
because `steps/` has no `__init__.py`, e.g. when a worktree's code is run from the main checkout's
folder). Only `run`, `verify`, `rejudge` and `papers-pilot` call the Claude API.
It reads `ANTHROPIC_API_KEY` from the environment or from a `.env` file in the current directory or the
repository root; `.env` is gitignored. Everything is resumable: the harvest saves its page token, taxonomy
answers and model answers are cached, and submitted batch ids are saved before the run waits on them.

## What each stage does

| Stage | What happens | Cost |
|---|---|---|
| `harvest` | Pages through the NCBI Datasets v2 genome reports for taxon 2 (Bacteria), GenBank, complete genome + chromosome, keeping the accession, organism, assembly facts, CheckM, NCBI's taxonomy check and every BioSample attribute (`ncbi/pages/*.jsonl.gz`). Then lists the metagenome-assembled genomes among them | free |
| `prepare` | Looks up every host name, host taxid and organism taxid in NCBI Taxonomy (about 30,000 lookups the first time, cached in `taxa_cache.json`). Looks up the ontology codes that decide a record's source on EBI OLS4 (about 400 terms, cached in `ontology_cache.json`; see "Ontology codes"). Then groups genomes into distinct BioSample records, applies the rules tier and estimates what the models would cost (`reports/prepare_report.md`) | free |
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

### Ontology codes

Some BioSamples give the sample only as an ontology code, alone or with its label: env_medium `ENVO:00002262`,
"clay soil [ENVO:00002262]". The vocabulary lists one id per term, so a code that is not a term's own id used to
go to the models. Now the rules place it through the ontology (`catalogue/ontology.py`). The code's terms are
looked up on EBI OLS4 during `prepare`, then:

1. If the code's own label or one of its synonyms is a vocabulary synonym, that is its term (ocean water,
   `ENVO:00002149`, is seawater).
2. Otherwise the code climbs to the nearest specific term within three steps (clay soil is a soil). ENVO and
   FoodOn codes climb is_a links only. In ENVO, a part_of link joins a material to where it lies, so estuarine
   mud is part of an estuary, not estuarine water. UBERON and PO codes also climb part_of links (mucosa of
   vagina is part of the vagina).

Climbing is held back in three cases, and the models then read the record as before:

- It never reaches a broad "type not stated" term. Estuarine mud is a kind of environmental material, but a
  model reading it says sediment.
- It never reaches a term listed under `not_placed_by_descent` in the vocabulary. ENVO's house includes cow
  sheds.
- It is not used when another source field of the record holds free text. In one record, env_medium codes
  "orchard" while plant_body_site says "leaf".

A code is also left to the models when its label is not the code's own (a wrong id), when the nearest terms
tie, or when the codes in one field disagree.

Measured on 2026-11-23, 214 distinct codes decide a record's source, and 60 of them lead to a term. The rule now
settles 82 records (167 genomes) that the models had answered:

| Compared with the earlier answer (curation, or the Fable judge's) | Records |
|---|---|
| Same term | 54 |
| Same niche | 8 |
| Earlier answer unknown | 8 |
| Different niche | 12 |

Of the 12 that change niche, about 3 are better (clam food product: shellfish rather than food), 2 or 3 are
worse ("food material" rather than fermented food), and the rest are debatable. The debatable ones are mostly
records whose medium is soil or seawater, where the models had preferred a host sample because a host is named.

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

When the rules can read the sample, they also take a host from:
- a field that is just a host name ("bovine");
- NCBI's metagenome names: "canine gut metagenome" is a dog, and "soil metagenome" has no host;
- a plant part or the rhizosphere, which give the host "plant, species not stated" (Fable's audit found these
  plant hosts missing in 310 genomes).

A record that names a living host but has no source field at all gets the source "host sample, site not
stated" (or "plant, part not stated"), with the host field as its evidence. A host field that names an
environment instead ("soil", "water") gives that environment as the source. Its source category (human-,
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

### Every dispute, judged twice (2026-11-26)

`rejudge --sample 100000` then had Opus 5.5 judge every settled dispute it can judge as a third party, for
$5.91. That was 2,379 records with 2,568 disputed parts. Opus decided as Fable's judge did on 84% of parts
(host 75%, source 86%). It disagreed in both directions: Fable's judge alone took Fable's answer 273 times, and
Opus alone took it 112 times.

`rejudge --model claude-sonnet-5-5 --curated-by claude-opus-5-5` gave the 431 disputes Opus had curated to Sonnet
5.5 instead, for $0.74. These records are the hardest: the two Sonnet runs had disagreed on them before Opus
settled them. Sonnet decided as Fable's judge did on only 51% of their 461 parts.

`assemble` lowers to `low` the confidence of every Fable-judge decision a second judge did not share, and lists
those parts in `judge_contested`: the source of 926 genomes and the host of 390.

## Rule audit (`audit-rules`)

The rules settle half the catalogue (about 54,000 genomes) and the cross-check never saw those records. A wrong
rule is wrong for every genome it touches. `audit-rules` has Fable label every record the rules settle on their
own and that holds some text (7,876 records, 53,659 genomes), with the same prompt, effort and cache as `verify`.
Each disagreement with a rule goes to the Fable judge, but the rule's answer always stands; the audit only
reports. `reports/audit_rules_*.json` gives the outcomes per part, weighted by records and by genomes, and
`*_patterns.tsv` groups the disputes by the rule's evidence (field and value), most genomes first, with Fable's
answers and the judge's verdicts. That is the list a wrong rule shows up in.

Audit run on 2026-11-26 ($31.79). It is weighted by genomes and judged by the Fable judge, so it is another
model's view, not ground truth:

| | Genomes the rules decide | The judge sided against the rule |
|---|---|---|
| Host | 50,371 | 206 (0.4%) |
| Source | 41,771 | 1,429 (3.4%) |

Most source disputes are rules that are vague, not wrong. The rules take the first informative field, and when
that field is broad, a later field often says more:

- env_medium "water" (390 genomes) is lake, sea or river water in another field.
- sample_type "tissue sample" (169) is a plant part.
- isolation_source "Plant" (103) is a root.
- "human", "patient" and "animal" (about 70) name a body site elsewhere in the record.

Some rules are plainly wrong:

- isolation_source "hospital" means a clinical isolate, not the hospital environment (31 of 31 records judged).
- Oyster and shrimp are hosts, not food, in aquaculture studies.
- The wild goat, *Capra aegagrus*, maps to "other mammal" rather than goat.

The rules also miss hosts the text states:

- "canine" in metagenome_source (691 genomes are dogs).
- A rhizosphere sample whose host can only be a plant (310 genomes).

Fixes applied on 2026-11-26 (vocabulary rules_revision 3):
1. **Refinement.** A first field that gives only a broad term ("water, type not stated", "host sample, site
   not stated", "plant, part not stated", "tissue, type not stated", "environmental sample", "food") is refined
   by a later field the rules can read, if that field names a more specific term of the same kind. Water is
   refined only by a kind of water: env_medium "water" with env_local_scale "lake water" is lake water.
2. **Host names in metagenome names.** NCBI's metagenome names carry the host ("canine gut metagenome").
3. **Plant hosts.** A plant part or the rhizosphere gives the host "plant, species not stated".
4. **Goat.** *Capra aegagrus* is a goat.
5. **Environments in host fields.** A host field that names an environment, on a record with no source
   field, gives the source.
6. **"hospital".** It is no longer a synonym of the hospital environment. Fable read it as a clinical isolate
   when the record names a patient or a disease, and as the building otherwise, so the models read each such
   record in context.

Checked against the audit's verdicts, which were saved before the change, the fixes improve the host of 1,026
genomes and the source of 320, and make none worse. 267 genomes with "hospital" went to the models. Oyster and
shrimp stay as they were: in a record they are food as often as hosts.

The missed sources the audit reports are mostly not misses. 3,483 of 3,633 genomes are records that name only
a host, and the catalogue already labels those "host sample, site not stated" or "plant, part not stated", as
Fable does; the audit compares with the rules alone.

## Papers (`papers-pilot`, `papers-fill`)

NCBI links many genomes' BioProjects to PubMed papers. The pilot measures what reading those papers would add,
before anyone pays to read all of them (`catalogue/papers.py`).

1. **Sample.** Genomes are drawn in random order. Each one whose BioProject links a paper adds one of its
   papers, until there are 200 distinct papers. Because the sample is over genomes, its gain rate is a rate
   over genomes.
2. **Texts.** The pilot uses the PubMed Central full text where NCBI serves it. Methods-like sections and the
   tables come first, then the other sections and the figure captions. Otherwise it uses the PubMed title and
   abstract.
3. **Read.** Claude Fable 5.2 reads each paper at high effort, one Batch API request per paper. It gets the
   text and up to 30 of the genomes NCBI links to the paper (accession, strain and isolate names); the
   sampled genomes are always among them. For the study as a whole, and for each genome the paper names
   individually, it reports eight facts: host, source, country, collection year, disease, isolation method,
   outbreak and phenotype.
   - Host and source must come from the catalogue vocabulary.
   - Every value needs a verbatim quote from the paper, and a value whose quote is not in the text is dropped.
4. **Compare.** Each genome's facts (its own where the paper names it, else the study's) are compared with
   its catalogue row:
   - `new`: the catalogue lacks the fact.
   - `same`, or `same niche` for a source.
   - `differs`: the paper and the catalogue disagree.

Pilot run on 2026-11-23 (seed 1; $34.52 in all, with a test request, the retries and the narrow request). Of
1,185 random genomes,
287 (24%) had a linked paper, which gave 200 distinct papers:

- **Texts:** 121 full texts, 77 abstracts, and 2 with neither.
- **Answers:** Fable 5.2 answered 158 of the 198 papers. It refused 39 and ran out of room on 1.
- **Refusals:** the refused papers are about pathogens and toxins (botulinum neurotoxin, superantigens, Shiga
  toxin, *Brucella*, *B. pseudomallei*, *B. anthracis*, cholera toxin phage, *M. tuberculosis*,
  *C. difficile*). `--retry-model claude-opus-5-5` got answers for only 4 of those 40.
- **Narrow request:** `--narrow-refused` then asked the remaining 36 for host, source, country and year only,
  from just the methods sections and tables. This cost $4.06, and 16 papers were answered: *Brucella* from 57
  patients, *B. pseudomallei*, multidrug-resistant *M. tuberculosis*, a Shiga toxin phage and streptococcal
  superantigens among them. Twenty papers (10% of the sample) are still refused. They include botulinum
  neurotoxin, *B. anthracis*, papers on genetic methods (knockout minicircles, resistant *Chlamydia*
  recombinants) and some ordinary population studies (*H. pylori*, *B. pertussis*). So the narrow request
  lifts the answered share from 82% to 90%.
- **Quotes:** every fact but one was found verbatim in its paper.

What the 223 sampled genomes on papers answered before the narrow request gained (the narrow request added 20
genomes; 3 of them gained a host, source, country or year):

| Gain over the catalogue | All | Full text | Abstract only |
|---|---|---|---|
| Any new fact | 74% | 83% | 62% |
| A host, source, country or year the catalogue lacked | 12% | 17% | 6% |
| Disease, isolation method, outbreak or phenotype | 71% | 80% | 60% |

By fact, the paper added a host for 18 of the 223 genomes, a source for 19, a country for 9 and a year for 7.
It added a phenotype for 90, an isolation method for 61, an outbreak or study context for 54, and a disease
for 15. Only 31% of these genomes are named individually in their paper; the rest get the paper's study-wide
facts. Those are often coarser than the genome's own BioSample: "bullfrog farms" against "feces".

Disagreements with the catalogue:

- **Source:** 42 genomes (19%). These are mixed: sometimes the paper is more specific (a bloodstream-infection
  cohort, where the catalogue says host sample), sometimes the BioSample is.
- **Host:** 5 genomes, including real BioSample errors. One paper's isolates are "clinical samples collected
  in 2021 from dogs", where BioSample says human; another's come from hospital sewage.
- **Year:** 9 genomes. Another 26 are a study period that contains the catalogue's year.

At these rates, reading every linked paper (about 9,850; 24% of genomes have one) would give about 2,500
genomes (2.3% of the catalogue) a host, source, country or year, and about 14,500 genomes free-text context.
At the pilot's $0.15 a paper, that is about $1,500 with Fable at high effort. Full texts carry most of the
core gain, for $0.19 a paper against $0.10 for an abstract.

`--dry-run` samples, fetches the texts and prices the requests without calling the API. The requests go out
in waves of 100. Each wave is checked against the budget at its worst case, every request using all of
`max_tokens`, so the budget holds even though a batch cannot be stopped part-way. Links, texts and batch ids
are cached in `papers/` and `papers_batches.json`. The report is `reports/papers_pilot_*.json`, and
`reports/papers_pilot_*.tsv` has one row per linked genome.

### Targeted run (`papers-fill`)

`papers-fill` reads the papers linked to every genome whose host or source the catalogue does not know, so that
`assemble` can fill those fields. It sends, budgets, caches and handles refusals as the pilot does, and keeps its
answers in `papers/answers_fill/`. Each request lists only the genomes to fill, up to 30 per paper. Options:

- `--full-text-only` skips papers that have only an abstract.
- `--core-only` asks only for host, source, country and year.
- `--narrow-refused` asks refused papers the narrow way.
- `--dry-run` prices the run for free.

Priced on 2026-11-26:

| | |
|---|---|
| Genomes with an unknown host or source | 17,167 |
| Of those, genomes with a linked paper | 4,673 (27%) |
| Distinct papers | 2,760: 1,129 full texts, 1,589 abstracts, 42 with neither |
| Expected cost, at the pilot's measured cost per paper ($0.222 a full text, $0.098 an abstract) | about $405 |
| Expected cost, full texts only | about $250 |
| Worst case, every request using all of its output allowance | $775 |

At the pilot's rates, reading these papers would give a host or source to roughly 1,500–2,200 more genomes.

### Merged into the catalogue

`assemble` reads every saved paper answer, from `papers-pilot` and `papers-fill`. For each genome a paper's
request listed, it fills the host, source, country and collection year that the catalogue does not know, and
never replaces a known value.

- **Several papers for one genome:** a paper that names the isolate wins over a study-wide statement, then full
  text over methods-only text over an abstract, then the lowest PMID.
- **How a filled value is marked:** a filled host or source gets tier `paper` (medium confidence) or
  `paper (study)` (low confidence). Its evidence field is `PMID <id>` and its quote is the paper's. Its
  category, niche and MIxS extension follow from the term, as for any other tier.
- **Disagreements:** these are flagged in `paper_conflicts` and never applied. Every fact also goes to
  `paper_facts.tsv.gz`.

With the pilot's 178 answered papers (2026-11-25), papers gave facts for 1,701 genomes. They filled 64
sources, 59 hosts, 35 countries and 25 collection years. The 491 host or source facts that disagree with the
catalogue were left as flags.

## Population structure (`structure-pilot`)

A yardstick for the ecotype methods (EcoSim, TreeCluster, fastbaps, PopPUNK) is a set of genomes whose
published group is known. `structure-pilot` builds one per species (`catalogue/structure.py`).

1. **Papers.** It searches PubMed for the species' population structure: reviews with the species in the title
   first, then other papers. A genus with only one named species in the catalogue is searched by the genus, so
   papers on "Prochlorococcus" count. It reads 15 papers: abstracts, plus the PubMed Central full text of two.
2. **Read.** In one Batch request per species, Claude Fable 5.2 lists the groups the papers define:
   subspecies, genomovars, lineages, phylogroups, clades, clonal complexes, sequence types, serotypes and
   ecotypes. For each it gives how the groups nest, what defines them, the reference strains named in them,
   and the ecology the papers report. Every group needs a verbatim quote from a cited paper.
3. **Match.** The reference strains are matched to the species' catalogue genomes by strain or isolate name.
   The species name, the type-strain mark and `=` aliases are handled: "KT2440 (mt-2 = DSM 6125)" is either
   name.

A species refused with full texts can be asked again from the abstracts alone with `--retry-refused`. The
report is `reports/structure_pilot_*.json`; `*_groups.tsv` holds the groups and `*_genomes.tsv` the matched
genomes.

Pilot run on 2026-11-27 for five species, $1.39 in all:

| Species | Read | Groups | Reference strains matched | Genomes with a known group |
|---|---|---|---|---|
| *Prochlorococcus marinus* | 15 papers, 2 in full | 7: high- and low-light ecotypes, with HLI, HLII, LLVII and LLVIII | 4 of 10 | 3 |
| *Pseudomonas putida* | 15 papers, 2 in full | 25: genomospecies, clonal complexes, sequence types | 20 of 69 | 26 |
| *Bacillus subtilis* | 15 papers, 2 in full | 5: subspecies *subtilis* and *spizizenii*, putative ecotype PE15 | 2 of 4 | 12 |
| *Escherichia coli* | 15 abstracts (the full texts were refused) | 22: phylogroups A–F, sequence types such as ST131 and ST1193, serogroups | none named | 0 |
| *Listeria monocytogenes* | refused, with full texts and from abstracts alone | | | |

Every quote was found in its paper. The groups and their ecology are sound, but few genomes get a known group,
because the papers name only a handful of reference strains each.

For typed pathogens, in-silico typing gives every genome its group without a model and without refusals:
MLST sequence types and clonal complexes, or ClermonTyping for *E. coli* phylogroups. The papers then say which
of those groups matter ecologically. Paper reading remains the only yardstick for species without a typing
scheme, such as the ecotypes of *Prochlorococcus* and *B. subtilis*.

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
| `*_evidence_field`, `*_quote`, `*_tier`, `*_model`, `*_confidence` | For host and source: the field and the exact words each value rests on, the tier that decided it, the model, and the tier's confidence. Tiers and confidences: `rule` high; `haiku`, `sonnet x2` and `fable judge` medium; `opus` low; `paper` medium and `paper (study)` low; `withheld`, `withheld by judge`, `not curated` and `no evidence` none. A paper's evidence field is `PMID <id>` |
| `host_crosscheck`, `source_crosscheck`, `fable_host`, `fable_source_term`, `crosscheck_model` | The cross-check outcome per part (see "Cross-check"; `rules only` where the rules decided the whole record, `not checked` before `verify`), Fable's own answers, and the Fable model and effort |
| `judge_contested` | The parts (`host`, `source`) whose Fable-judge decision a second judge did not share; their confidence is lowered to `low` (see "A second judge") |
| `paper_pmids`, `paper_filled`, `paper_conflicts` | The linked papers read for this genome (see "Papers"), which unknown fields they filled (`host`, `source`, `country`, `collection_year`), and which fields a paper disagrees with |
| `country`, `country_insdc`, `collection_year` | Parsed from `geo_loc_name` and `collection_date`; the country checked against the INSDC list. Filled from a paper only where BioSample gives none (see `paper_filled`) |
| `raw_*` | The BioSample values the curation read |
| `vocabulary_version` | Version of `steps/metadata_vocabulary.json` |

`paper_facts.tsv.gz`, one row per fact a paper gave for a genome:
- the PMID, the kind of text read, and the model;
- whether the paper names the isolate (`genome`) or speaks for all its isolates (`study`);
- the fact, its value and the quote;
- the catalogue's own value, and how the two compare: `new`, `same`, `same niche`, `within range`,
  `differs` or `in BioSample`;
- whether the fact filled the catalogue.

It also holds the facts the table has no column for (disease, isolation method, outbreak, phenotype),
and every disagreement, for review.

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
- 154 of the 214 ontology codes that decide a record's source are not placed, and those records keep the
  models' answers:
  - 116 reach no vocabulary term within three steps. Most are places and biomes such as ocean, stream,
    wetland area or agricultural field.
  - 31 reach only a broad or excluded term, such as rock, beach sand, estuarine mud or cow shed.
  - 7 are malformed or retired ids that OLS does not know.
- The rules take an environment field over a named host in 123 records (285 genomes), such as `env_medium`
  "seawater" on a marine sponge. Some of these are right (soil sampled around a named plant) and some are not.
- The vocabulary has no host group for protists or algae (ciliates, dinoflagellates, diatoms, amoebae). Their
  endosymbionts get an unknown host, or a wrong "animal" or "plant" one where the judge did not see them.
- For metagenome-assembled genomes the labels describe the metagenome sample, not an isolate.
- The vocabulary is fixed and coarse for some environments; free-text detail stays in the `raw_*` columns and
  the quotes.
