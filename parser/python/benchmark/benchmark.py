import os, sys
import time

p = os.path.abspath('..')
sys.path.append(p)

from gene_tree import GeneTree
import timeit
import random

from ete3 import PhyloxmlTree

from matplotlib import pyplot as plt


def run_listing_leaves(gt: GeneTree, n=10000):
    print("Test 2: Get all leaves")
    leaves = gt.trees[0].get_leaves()
    get_leaves_time = timeit.timeit(lambda: gt.trees[0].get_leaves(), number=n)
    print("Time ({} executions): {}".format(n, get_leaves_time))
    print("Time (Avg): {}".format(get_leaves_time / n))
    return leaves


def rand_homology_inference(gt: GeneTree, leaves, verbose=False):
    leaves = gt.trees[0].get_leaves()
    rand_pair = random.sample(leaves, 2)
    homology_type = gt.get_homology_type(rand_pair[0].name, rand_pair[1].name)
    if verbose:
        print("{}, {} {}".format(rand_pair[0].name, rand_pair[1].name, homology_type))
    return homology_type


def run_homology_inference(gt: GeneTree, leaves, n=10000, verbose=False):
    print("Test 3: Get homology type")
    homology_type_time = timeit.timeit(lambda: rand_homology_inference(gt, leaves, verbose=verbose), number=n)
    print("Time ({} executions): {}".format(n, homology_type_time))
    print("Time (Avg): {}".format(homology_type_time / n))
    print("Speed: {} ops/s".format(n / homology_type_time))
    return homology_type_time


g_temp = GeneTree()


def parse(file: str):
    g_temp.load_phylo_xml(file)


def run_parsing(file: str, n=1000):
    print("Test 1: Parsing")
    parsing_time = timeit.timeit(lambda: parse(file), number=n)
    file_size = os.path.getsize(file)
    gene_number = len(g_temp.trees[0].get_leaves())
    print("Time ({} executions): {}".format(n, parsing_time))
    print("Speed: {} MB/s".format(n * file_size / (parsing_time * 2**20)))
    print("Speed: {} genes/s".format(n * gene_number / parsing_time))


def run_finding_all_orthologs(gt: GeneTree, leaf: str, n=10000):
    print("Test 4: Getting all orthologs")
    getting_all_orthologs_time = timeit.timeit(lambda: gt.get_all_orthologs(leaf), number=n)
    gene_num = len(gt.trees[0].get_leaves())
    print("Tree size: {} genes".format(gene_num))
    print("Time ({} executions): {}".format(n, getting_all_orthologs_time))
    print("Time (Avg): {}".format(getting_all_orthologs_time / n))
    print("Speed: {} ops/s for gene tree of size {}".format(n / getting_all_orthologs_time, gene_num))


def run_finding_all_orthologs_indexed(gt: GeneTree, leaf: str, n=10000):
    print("Test 5: Getting all orthologs using interval index")
    gt.create_interval_index()
    getting_all_orthologs_time = timeit.timeit(lambda: gt.get_all_orthologs_indexed(leaf), number=n)
    gene_num = len(gt.trees[0].get_leaves())
    print("Tree size: {} genes".format(gene_num))
    print("Time ({} executions): {}".format(n, getting_all_orthologs_time))
    print("Time (Avg): {}".format(getting_all_orthologs_time / n))
    print("Speed: {} ops/s for gene tree of size {}".format(n / getting_all_orthologs_time, gene_num))


def plot_ortholog_inference_performance():
    gt1 = GeneTree()
    gt2 = GeneTree()
    gt3 = GeneTree()
    gt4 = GeneTree()
    gt5 = GeneTree()
    gt1.load_phylo_xml("test_data/BRCA_gene_tree_with_events.xml")
    gt2.load_phylo_xml("test_data/ENSACLG00000015576_with_events.xml")
    gt3.load_phylo_xml("test_data/ENSCSAVG00000010931_with_events.xml")
    gt4.load_phylo_xml("test_data/HOXA10_with_events.xml")
    gt5.load_phylo_xml("test_data/RNU6_92P_with_events.xml")
    size1, time1 = len(gt1.trees[0].get_leaves()), 100 / timeit.timeit(lambda: gt1.get_all_orthologs("ENSG00000139618"), number=100)
    size2, time2 = len(gt2.trees[0].get_leaves()), 100 / timeit.timeit(lambda: gt2.get_all_orthologs("ENSACLG00000015576"), number=100)
    size3, time3 = len(gt3.trees[0].get_leaves()), 100 / timeit.timeit(lambda: gt3.get_all_orthologs("ENSCSAVG00000010931"), number=100)
    size4, time4 = len(gt4.trees[0].get_leaves()), 100 / timeit.timeit(lambda: gt4.get_all_orthologs("ENSG00000253293"), number=100)
    size5, time5 = len(gt5.trees[0].get_leaves()), 100 / timeit.timeit(lambda: gt5.get_all_orthologs("ENSG00000272393"), number=100)
    plt.plot([size1, size2, size3, size4, size5], [time1, time2, time3, time4, time5], 'ro')
    plt.xlabel("Number of genes")
    plt.ylabel("ops/s")
    plt.title("Ortholog inference speed")
    plt.savefig("ortholog_inference_performance_number.png")


def plot_ortholog_inference_indexed_performance():
    gt1 = GeneTree()
    gt2 = GeneTree()
    gt3 = GeneTree()
    gt4 = GeneTree()
    gt5 = GeneTree()
    gt1.load_phylo_xml("test_data/BRCA_gene_tree_with_events.xml")
    gt2.load_phylo_xml("test_data/ENSACLG00000015576_with_events.xml")
    gt3.load_phylo_xml("test_data/ENSCSAVG00000010931_with_events.xml")
    gt4.load_phylo_xml("test_data/HOXA10_with_events.xml")
    gt5.load_phylo_xml("test_data/RNU6_92P_with_events.xml")
    gt1.create_interval_index()
    gt2.create_interval_index()
    gt3.create_interval_index()
    gt4.create_interval_index()
    gt5.create_interval_index()
    size1, time1 = len(gt1.trees[0].get_leaves()), 100 / timeit.timeit(lambda: gt1.get_all_orthologs_indexed("ENSG00000139618"), number=100)
    size2, time2 = len(gt2.trees[0].get_leaves()), 100 / timeit.timeit(lambda: gt2.get_all_orthologs_indexed("ENSACLG00000015576"), number=100)
    size3, time3 = len(gt3.trees[0].get_leaves()), 100 / timeit.timeit(lambda: gt3.get_all_orthologs_indexed("ENSCSAVG00000010931"), number=100)
    size4, time4 = len(gt4.trees[0].get_leaves()), 100 / timeit.timeit(lambda: gt4.get_all_orthologs_indexed("ENSG00000253293"), number=100)
    size5, time5 = len(gt5.trees[0].get_leaves()), 100 / timeit.timeit(lambda: gt5.get_all_orthologs_indexed("ENSG00000272393"), number=100)
    plt.plot([size1, size2, size3, size4, size5], [time1, time2, time3, time4, time5], 'ro')
    plt.xlabel("Number of genes")
    plt.ylabel("ops/s")
    plt.title("Ortholog inference speed (indexed)")
    plt.savefig("ortholog_inference_performance_number_indexed.png")


def plot_ortholog_inference_performance_height():
    gt1 = GeneTree()
    gt2 = GeneTree()
    gt3 = GeneTree()
    gt4 = GeneTree()
    gt5 = GeneTree()
    gt1.load_phylo_xml("test_data/BRCA_gene_tree_with_events.xml")
    gt2.load_phylo_xml("test_data/ENSACLG00000015576_with_events.xml")
    gt3.load_phylo_xml("test_data/ENSCSAVG00000010931_with_events.xml")
    gt4.load_phylo_xml("test_data/HOXA10_with_events.xml")
    gt5.load_phylo_xml("test_data/RNU6_92P_with_events.xml")
    size1, time1 = gt1.height(), 100 / timeit.timeit(lambda: gt1.get_all_orthologs("ENSG00000139618"), number=100)
    size2, time2 = gt2.height(), 100 / timeit.timeit(lambda: gt2.get_all_orthologs("ENSACLG00000015576"), number=100)
    size3, time3 = gt3.height(), 100 / timeit.timeit(lambda: gt3.get_all_orthologs("ENSCSAVG00000010931"), number=100)
    size4, time4 = gt4.height(), 100 / timeit.timeit(lambda: gt4.get_all_orthologs("ENSG00000253293"), number=100)
    size5, time5 = gt5.height(), 100 / timeit.timeit(lambda: gt5.get_all_orthologs("ENSG00000272393"), number=100)
    plt.plot([size1, size2, size3, size4, size5], [time1, time2, time3, time4, time5], 'ro')
    plt.xlabel("Gene tree height")
    plt.ylabel("ops/s")
    plt.title("Ortholog inference speed")
    plt.savefig("ortholog_inference_performance_height.png")


def plot_tree_height():
    gt1 = GeneTree()
    gt2 = GeneTree()
    gt3 = GeneTree()
    gt4 = GeneTree()
    gt5 = GeneTree()
    gt6 = GeneTree()
    gt1.load_phylo_xml("test_data/BRCA_gene_tree_with_events.xml")
    gt2.load_phylo_xml("test_data/ENSACLG00000015576_with_events.xml")
    gt3.load_phylo_xml("test_data/ENSCSAVG00000010931_with_events.xml")
    gt4.load_phylo_xml("test_data/HOXA10_with_events.xml")
    gt5.load_phylo_xml("test_data/RNU6_92P_with_events.xml")
    gt6.load_phylo_xml("test_data/HLA_F_with_events.xml")
    size1, height1 = len(gt1.trees[0].get_leaves()), gt1.height()
    size2, height2 = len(gt2.trees[0].get_leaves()), gt2.height()
    size3, height3 = len(gt3.trees[0].get_leaves()), gt3.height()
    size4, height4 = len(gt4.trees[0].get_leaves()), gt4.height()
    size5, height5 = len(gt5.trees[0].get_leaves()), gt5.height()
    size6, height6 = len(gt6.trees[0].get_leaves()), gt6.height()
    plt.plot([size1, size2, size3, size4, size5, size6], [height1, height2, height3, height4, height5, height6], 'bo')
    plt.xlabel("Number of genes")
    plt.ylabel("Tree height")
    plt.title("Tree height vs number of genes")
    plt.savefig("tree_height.png")


if __name__ == "__main__":
    gt = GeneTree()
    gt.load_phylo_xml("test_data/RNU6_92P_with_events.xml")
    run_parsing("test_data/RNU6_92P_with_events.xml", n=10)
    print("-" * 20)
    leaves = run_listing_leaves(gt)
    print("-" * 20)
    run_homology_inference(gt, leaves)
    print("-" * 20)
    run_finding_all_orthologs(gt, "ENSG00000272393", n=100)
    print("-" * 20)
    run_finding_all_orthologs_indexed(gt, "ENSG00000272393", n=100)
    plt.figure(0)
    plot_tree_height()
    plt.figure(1)
    plot_ortholog_inference_performance()
    plt.figure(2)
    plot_ortholog_inference_performance_height()
    plt.figure(3)
    plot_ortholog_inference_indexed_performance()
