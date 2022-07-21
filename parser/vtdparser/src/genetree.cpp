#pragma once

#include "genetree.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include <sstream>
#include <iostream>
#include <fstream>
#include "wchar.h"
#include "sys/stat.h"
#include <set>

using namespace compara;
using namespace vtdxml;
using namespace std;

GeneTree::GeneTree(const char* filename) {
    this->filename = filename;
    this->parse();
}

GeneTree::~GeneTree() {
    // TODO: destructor is buggy
    //delete this->vn;
    //delete this->root;
}

void GeneTree::parse() {
    FILE *f = NULL;
    unsigned char *xml = NULL;
    f = fopen(this->filename, "r");
    struct stat s;
    stat(this->filename, &s);
    xml = (unsigned char*)malloc(sizeof(unsigned char) * (int)s.st_size);
    fread(xml, sizeof(unsigned char), s.st_size, f);

    VTDGen *vg = NULL;
    VTDNav *vn = NULL;
    vg = new VTDGen();
    vg->setDoc_BR(xml, s.st_size);

    try {
        vg->parse(true);
        vn = vg->getNav();
        this->vn = vn;
        this->parse_genetree_node();
    } catch (VTDException& e) {
        printf("error message is %s \n", e.getMessage());
    }
    
}

void GeneTree::parse_genetree_node() {
    // begin from root
    this->vn->toNode(ROOT);
    // check for phyloxml root tag
    if (this->vn != NULL && this->vn->getCurrentIndex() != -1) {
        wstring rt_str = this->vn->toStringLowerCase(this->vn->getCurrentIndex());
        if (rt_str.find(L"phyloxml") == wstring::npos) {
            cerr << "Not a valid xml file." << endl;
            return;
        }
    }
    // traverse the immediate children to locate the phylogeny tag
    this->vn->toNode(FIRST_CHILD);
    int phylogeny_root = this->vn->getCurrentIndex();
    if (phylogeny_root != -1) {
        int curr = this->vn->getCurrentIndex();
        wstring curr_str = this->vn->toStringLowerCase(curr);
        while (curr != -1 && curr_str.find(L"phylogeny") == wstring::npos) {
            this->vn->toNode(NEXT_SIBLING);
            curr = this->vn->getCurrentIndex();
            curr_str = this->vn->toStringLowerCase(curr);
        }
        if (curr == -1) {
            cerr << "Did not find phylogeny tag." << endl;
            return;
        }
    } else {
        cerr << "Did not find root tag." << endl;
        return;
    }
    // traverse the children to locate the first clade tag
    this->vn->toNode(FIRST_CHILD);
    int clade_root = this->vn->getCurrentIndex();
    if (clade_root != -1) {
        int curr = this->vn->getCurrentIndex();
        wstring curr_str = this->vn->toStringLowerCase(curr);
        while (curr != -1 && curr_str.find(L"clade") == wstring::npos) {
            this->vn->toNode(NEXT_SIBLING);
            curr = this->vn->getCurrentIndex();
            curr_str = this->vn->toStringLowerCase(curr);
        }
        if (curr == -1) {
            cerr << "Did not find clade tag." << endl;
            return;
        }
    } else {
        cerr << "Did not find root tag." << endl;
        return;
    }
    this->root = new GeneTreeNode(new BookMark(this->vn->cloneNav()));
}

vector<wstring> GeneTree::get_genes() {
    vector<GeneTreeNode*> descendants = this->root->get_descendants();
    vector<wstring> genes;
    for (int i = 0; i < descendants.size(); i++) {
        GeneTreeNode *node = descendants[i];
        if (node->is_leaf()) {
            genes.push_back(node->get_name());
        }
    }
    return genes;
}

void GeneTree::print() {
    if (this->root != NULL) {
        this->root->print();
    }
}

/**
 * @brief Write the interval-based index of the gene tree to a file.
 * 
 * The index contains the following information: leaves and their labels,
 * internal nodes and their labels, and the labels of all internal nodes
 * that are duplication nodes.
 * 
 * @param filename 
 */
void GeneTree::write_index(const char *filename) {
    ofstream fout(filename, ios::out | ios::binary);
    vector<GeneTreeNode*> descendants = this->root->get_descendants();
    descendants.push_back(this->root);
    vector<GeneTreeNode*> genes;
    for (int i = 0; i < descendants.size(); i++) {
        GeneTreeNode *node = descendants[i];
        if (node->is_leaf()) {
            genes.push_back(node);
        }
    }
    int gene_num = genes.size();
    // Section 1: leaves and their labels
    // number of genes (leaves)
    fout.write((char*)&gene_num, sizeof(int));
    map<int, int> leaf_labels;
    for (int i = 0; i < genes.size(); i++) {
        genes[i]->write_index(i, fout);
        leaf_labels.insert(pair<int, int>(genes[i]->bm->hashCode(), i));
    }
    // Section 2: internal nodes and their labels
    // number of internal nodes
    int internal_nodes_num = descendants.size() - genes.size();
    int duplication_nodes_num = 0;
    vector<GeneTreeNode*> duplication_nodes;

    fout.write((char*)&internal_nodes_num, sizeof(int));
    for (int i = 0; i < descendants.size(); i++) {
        GeneTreeNode *node = descendants[i];
        if (!node->is_leaf()) {
            if (node->node_type == GeneTreeNodeType::DUPLICATION) {
                bool has_leaf = false;
                for (int j = 0; j < node->children.size(); j++) {
                    if (node->children[j]->is_leaf()) {
                        has_leaf = true;
                        break;
                    }
                }
                if (true) {
                    duplication_nodes.push_back(node);
                    duplication_nodes_num++;
                }
            }
            vector<GeneTreeNode*> leaves = node->get_leaves();
            int min_label = MAXINT;
            int max_label = MININT;
            // label of internal node is the min and max label of its leaves
            for (int j = 0; j < leaves.size(); j++) {
                GeneTreeNode *leaf = leaves[j];
                int curr_label = leaf_labels.at(leaf->bm->hashCode());
                if (curr_label < min_label) {
                    min_label = curr_label;
                }
                if (curr_label > max_label) {
                    max_label = curr_label;
                }
            }
            node->write_index(make_tuple(min_label, max_label), fout);
        }
    }

    // Section 3: duplication nodes and their labels
    // number of duplication nodes
    fout.write((char*)&duplication_nodes_num, sizeof(int));
    for (int i = 0; i < duplication_nodes.size(); i++) {
        GeneTreeNode *node = duplication_nodes[i];
        vector<GeneTreeNode*> leaves = node->get_leaves();
        int min_label = MAXINT;
        int max_label = MININT;
        // label of duplication node is the min and max label of its leaves
        for (int j = 0; j < leaves.size(); j++) {
            GeneTreeNode *leaf = leaves[j];
            int curr_label = leaf_labels.at(leaf->bm->hashCode());
            if (curr_label < min_label) {
                min_label = curr_label;
            }
            if (curr_label > max_label) {
                max_label = curr_label;
            }
        }
        node->write_index(make_tuple(min_label, max_label), fout);
    }
    fout.close();
}

void GeneTree::load_index(const char *filename) {
    ifstream fin(filename, ios::in | ios::binary);
    this->gti = GeneTreeIndex::read_genetree_index(fin);
    this->index_loaded = true;
    fin.close();
}

vector<OrthologPair> GeneTree::get_orthologs(string query_gene) {
    vector<OrthologPair> orthologs;
    vector<string> one_to_many_genes;
    if (this->index_loaded) {
        IndexedGeneTreeNode idx_gene_node = this->gti->leaves.at(query_gene);
        int node_hash = idx_gene_node.node_hash;
        GeneTreeNode *gene_node = this->root->leaves_map.at(node_hash);
        // gene tree nodes store more information so we use pointers
        // to prevent copying of the gene trees
        vector<GeneTreeNode*> ancestors = gene_node->get_ancestors();
        vector<IndexedGeneTreeNode> ancestors_idx;
        for (int i = 0; i < ancestors.size(); i++) {
            GeneTreeNode *ancestor = ancestors[i];
            int ancestor_hash = ancestor->bm->hashCode();
            IndexedGeneTreeNode ancestor_node = this->gti->internal_nodes.at(ancestor_hash);
            ancestors_idx.push_back(ancestor_node);
        }
        vector<int> visited;
        vector<int> one_to_one_labels;
        vector<int> one_to_many_labels;
        vector<int> many_to_many_labels;
        set<int> many_to_x_candidates;
        int duplication_on_path = 0;
        for (int j = 0; j < ancestors_idx.size(); j++) {
            int min_label = get<0>(ancestors_idx[j].internal_label);
            int max_label = get<1>(ancestors_idx[j].internal_label);
            if (ancestors_idx[j].node_type == SPECIATION) {
                vector<Interval<int,int>> subtree_dup_nodes = this->gti->find_duplication_subintervals(min_label, max_label);
                for (int k = 0; k < subtree_dup_nodes.size(); k++) {
                    Interval<int,int> dup_node = subtree_dup_nodes[k];
                    for (int l = max(dup_node.start, min_label); l <= min(dup_node.stop, max_label); l++) {
                        many_to_x_candidates.insert(l);
                    }
                }
            }
            if (ancestors_idx[j].node_type == DUPLICATION) {
                duplication_on_path++;
            }
            if (idx_gene_node.label >= min_label && idx_gene_node.label <= max_label) {
                // label is within the range
                for (int k = min_label; k <= max_label; k++) {
                    if (ancestors_idx[j].node_type == SPECIATION) {
                            // if the internal node is a speciation node, add all labels in the range that are not visited
                        if (k != idx_gene_node.label && std::find(visited.begin(), visited.end(), k) == visited.end()) {
                            if (many_to_x_candidates.find(k) != many_to_x_candidates.end()) {
                                if (duplication_on_path > 0) {
                                    many_to_many_labels.push_back(k);
                                } else {
                                    one_to_many_labels.push_back(k);
                                }
                            } else {
                                if (duplication_on_path > 0) {
                                    one_to_many_labels.push_back(k);
                                }
                                else {
                                    one_to_one_labels.push_back(k);
                                }
                            }
                        }
                    }
                    visited.push_back(k);
                }
            }
        }
        for (int i = 0; i < one_to_one_labels.size(); i++) {
            int label = one_to_one_labels[i];
            IndexedGeneTreeNode idx_node = this->gti->leaf_labels.at(label);
            OrthologPair ortho_pair = {
                .gene_name = query_gene,
                .taxon = "",
                .ortholog_name = idx_node.gene_name,
                .type = ONE_TO_ONE
            };
            orthologs.push_back(ortho_pair);
        }
        for (int i = 0; i < one_to_many_labels.size(); i++) {
            int label = one_to_many_labels[i];
            IndexedGeneTreeNode idx_node = this->gti->leaf_labels.at(label);
            OrthologPair ortho_pair = {
                .gene_name = query_gene,
                .taxon = "",
                .ortholog_name = idx_node.gene_name,
                .type = ONE_TO_MANY
            };
            orthologs.push_back(ortho_pair);
        }
        for (int i = 0; i < many_to_many_labels.size(); i++) {
            int label = many_to_many_labels[i];
            IndexedGeneTreeNode idx_node = this->gti->leaf_labels.at(label);
            OrthologPair ortho_pair = {
                .gene_name = query_gene,
                .taxon = "",
                .ortholog_name = idx_node.gene_name,
                .type = MANY_TO_MANY
            };
            orthologs.push_back(ortho_pair);
        }
    }
    return orthologs;
}