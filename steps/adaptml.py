"""AdaptML (Hunt, David, Gevers, Preheim, Alm & Polz 2008, Science 320:1081) re-implemented in Python 3.

A hidden Markov model on the tree: every node sits in one of K hidden "projected habitats", a habitat
changes along a branch of length t at one symmetric rate mu (stay: 1/K + (K-1)/K*exp(-K*mu*t)), and a
genome in habitat h is observed in environment e with probability E[h][e]. E and mu are learned, redundant
habitats are dropped, every node gets its most likely habitat, and clades that fit their habitat better
than label-shuffled data become populations.

    python -m steps.adaptml tree.nwk genome_metadata.csv out_dir [--outgroup outgroup] [--column niche]
                            [--clone-map clone_groups.csv] [--seed 42]

Inputs: the rerooted tree (rerooted_trees/<species>.nwk) and a CSV whose first column is the genome and
whose --column is its observed environment; niche_metadata.py's genome_metadata.csv fits as it is
(first column "accession", label column "niche" or "host").

Writes into out_dir:
  adaptml_membership.csv  genome -> population, in ecotype_membership.csv's columns (load_membership reads it).
                          Genomes outside every significant cluster have no row: AdaptML leaves them unassigned.
  adaptml_genomes.csv     every tree tip: environment, habitat, population (blank if unassigned)
  adaptml_habitats.csv    one row per learned habitat: its distribution over the environments
  adaptml_run.json        mu, log-likelihood, AIC and the parameters used

FIDELITY. The original is Python 2 at github.com/almlab/adaptml and has no license, so none of it is copied
here. The paper's supporting online material could not be downloaded (science.org refused), so the model was
taken from reading that source and its user guide (readme.pdf), NOT from the SOM. Nothing below has been
checked against the SOM text, and the original was never run, so outputs were not compared number for number.

Same as the original source:
  - transition model, random start of E, defaults (16 start habitats, collapse 0.10, converge 0.001, mu 1.0)
  - learning on the unrooted tree: sum-product partial likelihoods, flat prior, node posteriors
  - E update: each labelled tip adds, to every habitat h, the probability that the tip is in h predicted from
    its parent's posterior across the tip's branch (this is what the source does; it is not textbook EM)
  - mu update ("rateopt=avg"): share of internal branches whose two ends have different best habitats,
    inverted through the transition formula at the mean length of those branches
  - stop when E moves < converge_thresh per cell on average (or repeats, or after 502 rounds); then drop any
    habitat within squared distance collapse_thresh of one already kept and re-learn, until none is dropped
  - habitat of every node: joint maximum (Viterbi) on the tree rooted at the ingroup root, no root prior
  - a population is the topmost node of a same-habitat stretch where > 90% of the remaining tips have the
    node's habitat and the likelihood beats the cluster_thresh quantile of label-shuffled likelihoods (assigned
    habitats held fixed); deeper stretches are settled first and their clusters taken out of the ones above;
    a single tip is never a population; zero-length branches get the shortest positive length

Interpretation (the source is unclear, approximate, or crashes there):
  - after a nested cluster is taken out, the source patches the parent's likelihood and threshold by adding
    and subtracting per-node values; here both are recomputed over the tips that remain. Same test, exact null.
  - mu: the source divides by zero / takes log of <= 0 when no or all internal branches change habitat;
    here the value is clamped, and mu is kept when the tree has no internal branch
  - ties and the order habitats are compared in for collapsing: lowest index first (source: Python 2 dict order)
  - trees that are not binary are accepted (the model does not need binary; the source's parser does)

Different on purpose:
  - the outgroup is removed before learning. The source keeps it as a tip with an arbitrary environment label.
  - tips with no label (missing from the CSV, blank or "unknown") stay in the tree as missing data: they add
    nothing to E and take no part in shuffling, but still get a habitat. The source needs a label on every tip.
  - shuffles default to 1000 (source 10000); cluster_thresh default 0.95 is the source wrapper's default.
    Every node is tested at that level with no multiple-testing correction, as in the source, so random
    labels still yield a few small populations: read adaptml_habitats.csv before trusting them.

Left out: rateopt=num (numerical mu), the iTOL/colour outputs, tree truncation, cluster distances, and the lab
practice of 100 restarts to pick "stable habitats" (results depend on the random start: rerun with other seeds).
"""
import argparse
import csv
import json
import random
from math import exp, floor, inf, log
from pathlib import Path

from Bio import Phylo

from steps.parsing import load_clone_map

HOMOGENEITY = 0.9  # fixed in the original


def _flatten(root, unroot=False):
    """Preorder (tip names, parent index, child lists, branch lengths). unroot=True removes a two-child
    root by hanging one child off the other, which is the unrooted tree the habitats are learned on."""
    first = [(c, 0, c.branch_length or 0.0) for c in root.clades]
    if unroot and len(first) == 2:
        (a, _, la), (b, _, lb) = sorted(first, key=lambda item: item[0].is_terminal())
        first = [(c, 0, c.branch_length or 0.0) for c in a.clades] + [(b, 0, la + lb)]
    names, parent, length, stack = [None], [0], [0.0], first[::-1]
    while stack:
        clade, up, t = stack.pop()
        names.append(clade.name if clade.is_terminal() else None)
        parent.append(up)
        length.append(t)
        stack += [(c, len(names) - 1, c.branch_length or 0.0) for c in reversed(clade.clades)]
    shortest = min(t for t in length[1:] if t > 0)
    length = [t if t > 0 else shortest for t in length]
    kids = [[] for _ in names]
    for i in range(1, len(names)):
        kids[parent[i]].append(i)
    return names, parent, kids, length


def _across(v, stay):
    """Habitat vector v carried across a branch; stay = exp(-K*mu*t)."""
    base = (1 - stay) / len(v) * sum(v)
    return [base + stay * x for x in v]


def _scaled(v, by=max):
    top = by(v)
    return [x / top for x in v]


def _best(v):
    return max(range(len(v)), key=v.__getitem__)


def _round(tree, obs, E, mu):
    """One learning round on the unrooted tree: returns the updated (E, mu)."""
    parent, kids, length = tree
    k, n = len(E), len(parent)
    stay = [exp(-k * mu * t) for t in length]
    below, up = [None] * n, [None] * n  # evidence under a node; the same carried to its parent
    for i in reversed(range(n)):
        v = [1.0] * k
        if not kids[i] and obs[i] is not None:
            v = [row[obs[i]] for row in E]
        for c in kids[i]:
            v = [x * y for x, y in zip(v, up[c])]
        below[i] = _scaled(v)
        up[i] = _across(below[i], stay[i])
    above, post = [None] * n, [None] * n  # evidence from the rest of the tree; node posteriors
    above[0] = [1.0] * k
    for i in range(n):
        if not kids[i]:
            continue
        post[i] = _scaled([x * y for x, y in zip(above[i], below[i])], sum)
        for c in kids[i]:
            v = above[i]
            for sibling in kids[i]:
                if sibling != c:
                    v = [x * y for x, y in zip(v, up[sibling])]
            above[c] = _across(_scaled(v), stay[c])

    new = [[exp(-24)] * len(E[0]) for _ in E]
    for i in range(n):
        if not kids[i] and obs[i] is not None:
            for h, p in enumerate(_across(post[parent[i]], stay[i])):
                new[h][obs[i]] += p
    new = [_scaled(row, sum) for row in new]

    inner = [i for i in range(1, n) if kids[i]]  # branches with an internal node at both ends
    if inner:
        changed = sum(_best(post[i]) != _best(post[parent[i]]) for i in inner) / len(inner)
        left = min(max(1 - changed * k / (k - 1), 1e-9), 1 - 1e-9)  # clamp: the original crashes outside (0, 1)
        mu = -log(left) / (k * sum(length[i] for i in inner) / len(inner))
    return new, mu


def learn_habitats(tree, obs, n_env, init_habitats, collapse_thresh, converge_thresh, rng):
    """Returns (E, mu): E[h][e] for the habitats that survive collapsing."""
    E = [_scaled([rng.random() for _ in range(n_env)], sum) for _ in range(init_habitats)]
    mu = 1.0
    while True:
        diff = 1.0
        for _ in range(502):
            old_diff = diff
            new, mu = _round(tree, obs, E, mu)
            diff = sum(abs(a - b) for r1, r2 in zip(E, new) for a, b in zip(r1, r2)) / (len(E) * n_env)
            E = new
            if diff < converge_thresh or floor(diff * 1e8) == floor(old_diff * 1e8) > 0:
                break
        kept = []
        for row in E:
            if all(sum((a - b) ** 2 for a, b in zip(row, other)) >= collapse_thresh for other in kept):
                kept.append(row)
        if len(kept) == len(E) or len(kept) < 2:
            return kept, mu
        E = kept


def assign_habitats(tree, obs, E, mu):
    """Most likely joint habitat of every node of the rooted tree. Returns (habitat per node, log-likelihood)."""
    parent, kids, length = tree
    k, n = len(E), len(parent)
    score, choice = [None] * n, [None] * n
    for i in reversed(range(n)):
        if kids[i]:
            score[i] = [sum(score[c][h] for c in kids[i]) for h in range(k)]
        else:
            score[i] = [log(row[obs[i]]) if obs[i] is not None else 0.0 for row in E]
        if i:
            stay = exp(-k * mu * length[i])
            move = log((1 - stay) / k) if stay < 1 and k > 1 else -inf
            options = [[score[i][h2] + (log((1 - stay) / k + stay) if h == h2 else move) for h2 in range(k)]
                       for h in range(k)]
            choice[i] = [_best(row) for row in options]
            score[i] = [max(row) for row in options]
    habitat = [_best(score[0])] * n
    for i in range(1, n):
        habitat[i] = choice[i][habitat[parent[i]]]
    return habitat, max(score[0])


def find_populations(tree, obs, E, habitat, shuffles, cluster_thresh, rng):
    """Lists of tip indices, one per population, in tree order."""
    parent, kids, _ = tree
    n = len(parent)
    log_e = [[log(p) for p in row] for row in E]
    labelled = [i for i in range(n) if not kids[i] and obs[i] is not None]
    pool = [obs[i] for i in labelled]
    zeros = [0.0] * shuffles
    null = {i: [] for i in labelled}  # per tip: its log-likelihood in each shuffle, habitats held fixed
    for _ in range(shuffles):
        rng.shuffle(pool)
        for i, e in zip(labelled, pool):
            null[i].append(log_e[habitat[i]][e])

    def combine(parts):
        return (sum((p[0] for p in parts), []), sum(p[1] for p in parts),
                [sum(x) for x in zip(*(p[2] for p in parts))])

    # Per node (tips, log-likelihood, shuffled log-likelihoods) over the tips not already in a cluster:
    # test[i] when node i is tested, rest[i] once its own same-habitat stretch has been settled too.
    test, rest, members = [None] * n, [None] * n, {}
    for i in reversed(range(n)):
        if not kids[i]:
            test[i] = rest[i] = ([i], log_e[habitat[i]][obs[i]] if i in null else 0.0, null.get(i, zeros))
            continue
        tips, lik, shuffled = test[i] = combine([test[c] if habitat[c] == habitat[i] else rest[c] for c in kids[i]])
        if (tips and sum(habitat[t] == habitat[i] for t in tips) / len(tips) > HOMOGENEITY
                and lik > sorted(shuffled)[int(shuffles * cluster_thresh)]):
            members[i] = tips
            rest[i] = ([], 0.0, zeros)
        else:
            rest[i] = combine([rest[c] for c in kids[i]])
        for c in kids[i]:
            test[c] = rest[c] = None

    populations, covered = [], [False] * n  # covered: inside a cluster rooted higher in the same stretch
    for i in range(n):
        covered[i] = i > 0 and habitat[i] == habitat[parent[i]] and covered[parent[i]]
        if i in members and not covered[i]:
            populations.append(members[i])
            covered[i] = True
    return populations


def run_adaptml(tree_path, habitat_csv, out_dir, outgroup=None, column="niche", clone_map=None,
                init_habitats=16, collapse_thresh=0.10, converge_thresh=0.001,
                shuffles=1000, cluster_thresh=0.95, seed=42):
    """Run AdaptML on a rooted tree. habitat_csv: first column genome, `column` its environment.
    clone_map (dedup_clones.py's clone_groups.csv) gives collapsed genomes their representative's
    population. Returns the adaptml_membership.csv path."""
    with open(habitat_csv, newline="") as f:
        rows = list(csv.DictReader(f))
    if not rows or column not in rows[0]:
        raise ValueError(f"{habitat_csv} has no '{column}' column")
    key = next(iter(rows[0]))
    label = {r[key]: r[column].strip() for r in rows if r[column].strip().lower() not in ("", "unknown")}

    tree = Phylo.read(tree_path, "newick")
    if outgroup:
        tree.prune(outgroup)
    names, *rooted = _flatten(tree.root)
    envs = sorted({label[t] for t in names if t in label})
    if len(envs) < 2:
        raise ValueError(f"AdaptML needs at least two environments among the tree's genomes; found {envs}")
    obs = [envs.index(label[t]) if t in label else None for t in names]
    print(f"AdaptML: {sum(t is not None for t in names)} genomes, {sum(o is not None for o in obs)} with a "
          f"'{column}' label, {len(envs)} environments")

    rng = random.Random(seed)
    unrooted_names, *unrooted = _flatten(tree.root, unroot=True)
    E, mu = learn_habitats(unrooted, [envs.index(label[t]) if t in label else None for t in unrooted_names],
                           len(envs), init_habitats, collapse_thresh, converge_thresh, rng)
    habitat, loglik = assign_habitats(rooted, obs, E, mu)
    populations = find_populations(rooted, obs, E, habitat, shuffles, cluster_thresh, rng)
    population_of = {tip: number for number, tips in enumerate(populations, 1) for tip in tips}
    print(f"AdaptML: {len(E)} habitat(s), {len(populations)} population(s), mu = {mu:.4g}")

    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    tips = [i for i, t in enumerate(names) if t is not None]
    with open(out_dir / "adaptml_habitats.csv", "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["habitat", "genomes"] + envs)
        for h, row in enumerate(E):
            writer.writerow([f"H{h + 1}", sum(habitat[i] == h for i in tips)] + [f"{p:.4f}" for p in row])
    with open(out_dir / "adaptml_genomes.csv", "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["taxon", "environment", "habitat", "population"])
        for i in tips:
            writer.writerow([names[i], label.get(names[i], ""), f"H{habitat[i] + 1}", population_of.get(i, "")])
    (out_dir / "adaptml_run.json").write_text(json.dumps({
        "habitats": len(E), "populations": len(populations), "mu": mu, "log_likelihood": loglik,
        "aic": 2 * (1 + len(E) * (len(envs) - 1)) - 2 * loglik, "column": column, "seed": seed,
        "init_habitats": init_habitats, "collapse_thresh": collapse_thresh, "converge_thresh": converge_thresh,
        "shuffles": shuffles, "cluster_thresh": cluster_thresh}, indent=2) + "\n")

    clones = load_clone_map(clone_map)
    out = out_dir / "adaptml_membership.csv"
    with open(out, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["file", "gene_count", "trial", "ecotype", "ecotype_size", "taxon", "representative"])
        for number, members in enumerate(populations, 1):
            genomes = [(taxon, names[i]) for i in members for taxon in clones.get(names[i], [names[i]])]
            # gene_count left blank: load_membership int()s it when the file is not the EcoSim full-core XML
            writer.writerows(["adaptml", "", "", number, len(genomes), taxon, rep] for taxon, rep in genomes)
    print(f"Saved {out}")
    return out


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("tree")
    parser.add_argument("habitat_csv")
    parser.add_argument("out_dir")
    parser.add_argument("--outgroup")
    parser.add_argument("--column", default="niche")
    parser.add_argument("--clone-map")
    parser.add_argument("--seed", type=int, default=42)
    args = parser.parse_args()
    run_adaptml(args.tree, args.habitat_csv, args.out_dir, args.outgroup or None, args.column, args.clone_map,
                seed=args.seed)
