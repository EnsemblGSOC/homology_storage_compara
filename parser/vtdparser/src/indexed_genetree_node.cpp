#include "indexed_genetree_node.h"
#include "stdio.h"
#include "stdlib.h"
#include <string>
#include <sstream>
#include <iostream>
#include "wchar.h"

using namespace compara;
using namespace vtdxml;
using namespace std;

IndexedGeneTreeNode::IndexedGeneTreeNode(int node_hash,
                                         tuple<int, int> internal_label,
                                         GeneTreeNodeType node_type) {
    this->node_hash = node_hash;
    this->label = -1;
    this->internal_label = internal_label;
    this->gene_name = "";
    this->node_type = node_type;
}

IndexedGeneTreeNode::IndexedGeneTreeNode(int node_hash,
                                         int label,
                                         string gene_name) {
    this->node_hash = node_hash;                                        
    this->label = label;
    this->internal_label = make_tuple(-1, -1);
    this->gene_name = gene_name;
    this->node_type = GeneTreeNodeType::LEAF;
}

IndexedGeneTreeNode::~IndexedGeneTreeNode() {
}

string IndexedGeneTreeNode::get_name() {
    return this->gene_name;
}

void IndexedGeneTreeNode::write_index(ostream &out) {
    int hash = this->node_hash;
    // first 4 bytes: type of node represented by this line
    out.write((char*)&this->node_type, sizeof(int));
    if (this->node_type == GeneTreeNodeType::LEAF) {
        // second 4 bytes: label
        out.write((char*)&this->label, sizeof(int));
        // third 4 bytes: gene name length
        int gene_name_length = this->gene_name.length();
        out.write((char*)&gene_name_length, sizeof(int));
        // gene name (variable length)
        out.write(this->gene_name.c_str(), gene_name_length * sizeof(char));
    } else {
        // second 8 bytes: internal label
        out.write((char*)&get<0>(this->internal_label), sizeof(int));
        out.write((char*)&get<1>(this->internal_label), sizeof(int));
    }
    out.write((char*)&hash, sizeof(int));
}

