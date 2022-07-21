#pragma once

#include "vtdGen.h"
#include "vtdNav.h"
#include "vtdNav_L5.h"
#include "vtdException.h"
#include "bookMark.h"
#include <iostream>
#include <vector>
#include <tuple>
#include <map>

using namespace std;

namespace compara {
    /**
     * @brief Types of gene tree node
     */
    typedef enum {
        SPECIATION,
        DUPLICATION,
        DUBIOUS,
        GENE_SPLIT,
        LEAF,
        OTHER
    } GeneTreeNodeType;

    /**
     * @brief A gene tree node.
     */
    class GeneTreeNode {
        friend class GeneTree;
        public:
            GeneTreeNode(vtdxml::BookMark *bm);
            ~GeneTreeNode();
            wstring get_name();
            wstring get_element_type();
            GeneTreeNode* get_parent();
            vector<GeneTreeNode*> get_children();
            vector<GeneTreeNode*> get_ancestors();
            vector<GeneTreeNode*> get_descendants();
            vector<GeneTreeNode*> get_leaves();
            bool is_leaf();
            bool is_root();
            vector<GeneTreeNode*> children;
            // vector<GeneTreeNode*> ancestors;
            GeneTreeNode* parent;
            GeneTreeNodeType node_type;
            // A hashmap mapping from the hash of each node to the corresponding GeneTreeNode object
            map<int, GeneTreeNode*> leaves_map;
            void print(int depth = 0);
            void write_index(int label, ostream &out);
            void write_index(tuple<int, int>, ostream &out);
        
        private:
            void load_children();
            void load_node_type();
            wstring get_first_node(wstring node_name, wstring attrib_name = L"");

        protected:
            // default constructor, do not use
            GeneTreeNode();
            vtdxml::BookMark *bm;
    };
}
