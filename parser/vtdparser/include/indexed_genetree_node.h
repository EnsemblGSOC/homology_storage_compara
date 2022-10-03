#pragma once

#include "vtdGen.h"
#include "vtdNav.h"
#include "vtdNav_L5.h"
#include "vtdException.h"
#include "bookMark.h"
#include <iostream>
#include <vector>
#include <tuple>
#include "genetree_node.h"

using namespace compara;
using namespace std;

namespace compara {
    class IndexedGeneTreeNode {
        friend class GeneTree;
        friend class GeneTreeIndex;
        public:
            IndexedGeneTreeNode(int node_hash,
                                int label = -1,
                                string gene_name = "");
            IndexedGeneTreeNode(int node_hash,
                                tuple<int, int> internal_label,
                                GeneTreeNodeType node_type);
            ~IndexedGeneTreeNode();
            string get_name();
            GeneTreeNodeType node_type;

            void write_index(ostream &out);
            void read_index(istream &in);
            
        private:
            int node_hash;
            int label;
            tuple<int, int> internal_label;
            string gene_name;
    };
}