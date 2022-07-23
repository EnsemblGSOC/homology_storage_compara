#include "genetree_index.h"
#include "interval_tree.h"

using namespace std;
using namespace compara;

GeneTreeIndex::GeneTreeIndex() {
}

GeneTreeIndex* GeneTreeIndex::read_genetree_index(istream &in) {
    GeneTreeIndex *index = new GeneTreeIndex();
    index->load_leaves(in);
    index->load_internal_nodes(in);
    index->load_duplication_nodes(in);
    return index;
}

void GeneTreeIndex::load_leaves(istream &in) {
    int gene_num;
    int label;
    GeneTreeNodeType node_type;
    int gene_name_length;
    int node_hash;
    
    in.read((char*)&gene_num, sizeof(int));
    while (gene_num > 0 && !in.eof()) {
        in.read((char*)&node_type, sizeof(GeneTreeNodeType));
        if (node_type == GeneTreeNodeType::LEAF) {
            in.read((char*)&label, sizeof(int));
            in.read((char*)&gene_name_length, sizeof(int));
            char *gene_name = new char[gene_name_length + 1];
            in.read(gene_name, gene_name_length);
            gene_name[gene_name_length] = '\0';
            string gene_name_string = string(gene_name);
            in.read((char*)&node_hash, sizeof(int));
            IndexedGeneTreeNode node = IndexedGeneTreeNode(node_hash, label, gene_name_string);
            this->leaves.insert(pair<string, IndexedGeneTreeNode>(gene_name_string, node));
            this->leaf_labels.insert(pair<int, IndexedGeneTreeNode>(label, node));
            delete[] gene_name;
        }
        gene_num--;
    }
}

void GeneTreeIndex::load_internal_nodes(istream &in) {
    int internal_node_num;
    tuple<int, int> label;
    GeneTreeNodeType node_type;
    int node_hash;

    in.read((char*)&internal_node_num, sizeof(int));
    while (internal_node_num > 0 && !in.eof()) {
        in.read((char*)&node_type, sizeof(int));
        if (node_type != GeneTreeNodeType::LEAF) {
            int label_1;
            int label_2;
            in.read((char*)&label_1, sizeof(int));
            in.read((char*)&label_2, sizeof(int));
            in.read((char*)&node_hash, sizeof(int));
            label = make_tuple(label_1, label_2);
            this->internal_nodes.insert(pair<int, IndexedGeneTreeNode>(node_hash, IndexedGeneTreeNode(node_hash, label, node_type)));
        }
        internal_node_num--;
    }
}

void GeneTreeIndex::load_duplication_nodes(istream &in) {
    int duplication_node_num;
    int node_hash;
    int start;
    int end;
    GeneTreeNodeType node_type;
    IntervalTree<int, int>::interval_vector intervals;
    in.read((char*)&duplication_node_num, sizeof(int));
    while (duplication_node_num > 0 && !in.eof()) {
        in.read((char*)&node_type, sizeof(int));
        in.read((char*)&start, sizeof(int));
        in.read((char*)&end, sizeof(int));
        in.read((char*)&node_hash, sizeof(int));
        if (node_type == GeneTreeNodeType::DUPLICATION) {
            intervals.push_back(Interval<int, int>(start, end, node_hash));
        }
        duplication_node_num--;
    }
    this->duplication_nodes = IntervalTree<int, int>(move(intervals));
}

vector<Interval<int, int>> GeneTreeIndex::find_duplication_subintervals(int start, int end) {
    vector<Interval<int, int>> subintervals;
    IntervalTree<int, int>::interval_vector intervals = this->duplication_nodes.findContained(start, end);
    for (Interval<int, int> interval : intervals) {
        subintervals.push_back(interval);
    }
    return subintervals;
}

vector<IndexedGeneTreeNode> GeneTreeIndex::find_duplication_subtree_nodes(int node_hash) {
    vector<IndexedGeneTreeNode> subtree_nodes;
    IndexedGeneTreeNode query_node = this->internal_nodes.at(node_hash);
    int start = get<0>(query_node.internal_label);
    int end = get<1>(query_node.internal_label);
    if (start != -1 && end != -1) {
        vector<Interval<int, int>> subintervals = this->find_duplication_subintervals(start, end);
        for (Interval<int, int> interval : subintervals) {
            subtree_nodes.push_back(this->internal_nodes.at(interval.value));
        }
    }
    return subtree_nodes;
}