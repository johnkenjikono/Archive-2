"""Run Panaroo with a fix for its KeyError in collapse_families (panaroo <= 1.8.0).

merge_node_cluster de-duplicates the 'centroid' and 'protein' lists independently, but
downstream code zips them. When two centroids share an identical protein the lists
misalign, zip() drops a centroid, and centroid_to_index[...] raises KeyError.
Usage: python panaroo_patched.py <panaroo args>   (run inside the panaroo conda env)
"""
import sys
import panaroo.merge_nodes, panaroo.merge_graphs, panaroo.clean_network, panaroo.__main__

_orig = panaroo.merge_nodes.merge_node_cluster


def merge_node_cluster(G, nodes, newNode, *args, **kwargs):
    pairs = {}
    for n in nodes:
        for c, p in zip(G.nodes[n]['centroid'], G.nodes[n]['protein']):
            pairs.setdefault(c, p)
    G = _orig(G, nodes, newNode, *args, **kwargs)
    G.nodes[newNode]['protein'] = [pairs[c] for c in G.nodes[newNode]['centroid']]
    return G


# ponytail: monkeypatch every module that imported the name; drop this file once upstream fixes it
for mod in (panaroo.merge_nodes, panaroo.merge_graphs, panaroo.clean_network, panaroo.__main__):
    if hasattr(mod, 'merge_node_cluster'):
        mod.merge_node_cluster = merge_node_cluster

if __name__ == '__main__':
    sys.argv[0] = 'panaroo'
    sys.exit(panaroo.__main__.main())
