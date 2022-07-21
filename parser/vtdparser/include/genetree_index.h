#pragma once

#include "indexed_genetree_node.h"
#include "stdio.h"
#include "stdlib.h"
#include <string>
#include <sstream>
#include <iostream>
#include <map>
#include <tuple>
#include "wchar.h"
#include "interval_tree.h"

namespace compara {
    class GeneTreeIndex {
        public:
            static GeneTreeIndex *read_genetree_index(istream &in);
            // Map from gene name to IndexedGeneTreeNode containing node label
            map<string, IndexedGeneTreeNode> leaves;
            // Map from leaf label to IndexedGeneTreeNode with the label
            map<int, IndexedGeneTreeNode> leaf_labels;
            // Map from internal node hash to IndexedGeneTreeNode containing node label
            map<int, IndexedGeneTreeNode> internal_nodes;
            // Interval tree for duplication nodes
            IntervalTree<int, int> duplication_nodes;
            vector<IndexedGeneTreeNode> find_duplication_subtree_nodes(int node_hash);
        
        private:
            GeneTreeIndex();
            ~GeneTreeIndex();
            void load_leaves(istream &in);
            void load_internal_nodes(istream &in);
            void load_duplication_nodes(istream &in);
            vector<Interval<int, int>> find_duplication_subintervals(int start, int end);
    };
}