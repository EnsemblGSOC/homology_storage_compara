import glob
import os, sys
import timeit
from matplotlib import pyplot as plt
import random
p = os.path.abspath('..')
sys.path.append(p)
import gene_tree

paths = glob.glob('data_in/*.xml')
for p in paths:
    print("Processing {}".format(p))
    gt = gene_tree.GeneTree()
    gt.load_phylo_xml(p)
    gt.annotate_event_nodes_fast()
    gt.remove_seq_data()
    out_path = p.replace('data_in', 'data_out_without_seq')
    gt.export_phylo_xml(out_path)

out_paths = glob.glob('data_out_without_seq/*.xml')
sizes, speeds, speeds_indexed = [], [], []
for p in out_paths:
    print("Generating metrics for {}".format(p))
    gt = gene_tree.GeneTree()
    gt.load_phylo_xml(p)
    leaf_count = len(gt.trees[0].get_leaves())
    if leaf_count < 200:
        print("Small tree, skipping")
        continue
    gt.create_interval_index()
    sizes.append(leaf_count)
    leaves = [x.name for x in gt.trees[0].get_leaves()]
    time_unindexed = timeit.timeit(lambda: gt.get_all_orthologs(random.choice(leaves)), number=2000)
    time_indexed = timeit.timeit(lambda: gt.get_all_orthologs_indexed(random.choice(leaves)), number=2000)
    speeds.append(2000 / time_unindexed)
    speeds_indexed.append(2000 / time_indexed)
    print("{} seconds per operation indexed".format(time_indexed / 2000))
    print("{} genes, {} ops/s, {} ops/s indexed".format(sizes[-1], speeds[-1], speeds_indexed[-1]))
    print("{} times faster indexed".format(speeds_indexed[-1] / speeds[-1]))
    print("")

print(sizes)
print(speeds)
plt.plot(sizes, speeds, 'ro', label='Non-indexed')
plt.plot(sizes, speeds_indexed, 'bo', label='Indexed')
plt.xlabel("Number of genes")
plt.ylabel("ops/s")
plt.title("Ortholog inference speed (indexed)")
plt.legend(loc='upper right')
plt.savefig("ortholog_inference_performance_number_indexed.png")
