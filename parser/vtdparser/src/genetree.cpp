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
#include <unordered_set>
#include <unordered_map>

using namespace compara;
using namespace vtdxml;
using namespace std;

const tuple<int, int> DEFAULT_TUPLE = make_tuple(-1, -1);

GeneTree::GeneTree(const char* filename, bool is_vtdxml) {
    this->filename = filename;
    this->parse(is_vtdxml);
}

GeneTree::~GeneTree() {
    // TODO: destructor is buggy
    //delete this->vn;
    if (this->vn != NULL) {
        delete this->vn;
    }
    if (this->root != NULL) {
        delete this->root;
    }
}

void GeneTree::parse(bool is_vtdxml) {
    // FILE *f = NULL;
    // unsigned char *xml = NULL;
    // f = fopen(this->filename, "r");
    // struct stat s;
    // stat(this->filename, &s);
    // xml = (unsigned char*)malloc(sizeof(unsigned char) * (int)s.st_size);
    // fread(xml, sizeof(unsigned char), s.st_size, f);

    VTDGen *vg = NULL;
    VTDNav *vn = NULL;
    vg = new VTDGen();
    //vg->setDoc_BR(xml, s.st_size);
    if (!is_vtdxml) {
        try {
            vg->parseFile(true, this->filename);
            this->vn = vg->getNav();
            this->parse_genetree_node();
        } catch (...) {
            printf("Error: %s\n", this->filename);
            exit(1);
        }
    } else {
        try {
            this->vn = vg->loadIndex(this->filename);
            this->parse_genetree_node();
        } catch (...) {
            printf("Error: %s\n", this->filename);
            exit(1);
        }
    }
    
}

void GeneTree::write_vtdxml(const char* filename) {
    this->vn->writeIndex(const_cast<char*>(filename));
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
    unordered_map<int, int> leaf_labels;
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
                duplication_nodes.push_back(node);
                duplication_nodes_num++;
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
        unordered_set<int> visited;
        unordered_set<int> one_to_one_labels;
        unordered_set<int> one_to_many_labels;
        unordered_set<int> many_to_many_labels;
        int duplication_on_path = 0;
        
        vector<tuple<tuple<int, int>, tuple<int, int>>> path_labels;
        tuple<int, int> curr = ancestors_idx[0].internal_label;
        tuple<int, int> prev = make_tuple(idx_gene_node.label, idx_gene_node.label);

        for (int j = 0; j < ancestors_idx.size() - 1; j++) {
            if (ancestors_idx[j].node_type == SPECIATION || ancestors_idx[j].node_type == DUBIOUS) {
                if (curr != DEFAULT_TUPLE && prev != DEFAULT_TUPLE) {
                    path_labels.push_back(make_tuple(curr, prev));
                }
            } else if (ancestors_idx[j].node_type == DUPLICATION) {
                // include dummy duplication label for counting the number of duplications on the path
                path_labels.push_back(make_tuple(DEFAULT_TUPLE, DEFAULT_TUPLE));
            }
            prev = curr;
            curr = ancestors_idx[j+1].internal_label;
        }
        if (curr != DEFAULT_TUPLE && prev != DEFAULT_TUPLE) {
            if (ancestors_idx[ancestors_idx.size() - 1].node_type == SPECIATION || ancestors_idx[ancestors_idx.size() - 1].node_type == DUBIOUS) {
                path_labels.push_back(make_tuple(curr, prev));
            } else if (ancestors_idx[ancestors_idx.size() - 1].node_type == DUPLICATION) {
                path_labels.push_back(make_tuple(DEFAULT_TUPLE, DEFAULT_TUPLE));
            }
        }

        for (int k = 0; k < path_labels.size(); k++) {
            tuple<int, int> exclude = get<1>(path_labels[k]);
            tuple<int, int> include = get<0>(path_labels[k]);
            if (exclude != DEFAULT_TUPLE && include != DEFAULT_TUPLE) {
                // rightward path
                if (get<0>(exclude) == get<0>(include)) {
                    // get the set of duplication intervals in the right subtrees
                    vector<Interval<int, int>> sub_dup_nodes = this->gti->duplication_nodes.findContained(get<1>(exclude) + 1, get<1>(include));
                    int prev_min = 0, prev_max = 0;
                    for (int i = 0; i < sub_dup_nodes.size(); i++) {
                        int start = sub_dup_nodes[i].start;
                        int stop = sub_dup_nodes[i].stop;
                        // if the current interval is completely contained in a previous interval, skip it
                        if (start > prev_min && stop < prev_max) {
                            continue;
                        }
                        // trim the interval to remove the portion overlapping with a previous interval
                        if (start < prev_min && stop < prev_max && stop >= prev_min) {
                            stop = prev_min;
                        } else if (start > prev_min && stop > prev_max && start <= prev_max) {
                            start = prev_max;
                        }
                        for (int j = start; j <= stop; j++) {
                            if (duplication_on_path > 0) {
                                many_to_many_labels.insert(j);
                            } else {
                                one_to_many_labels.insert(j);
                            }
                            visited.insert(j);
                        }
                        prev_min = start;
                        prev_max = stop;
                    }

                    for (int i = get<1>(exclude) + 1; i <= get<1>(include); i++) {
                        if (visited.find(i) == visited.end()) {
                            if (duplication_on_path > 0) {
                                one_to_many_labels.insert(i);
                            } else {
                                one_to_one_labels.insert(i);
                            }
                        }
                    }
                }
                 // leftward path
                else if (get<1>(exclude) == get<1>(include)) {
                    // get the set of duplication intervals in the left subtrees
                    vector<Interval<int, int>> sub_dup_nodes = this->gti->duplication_nodes.findContained(get<0>(include), get<0>(exclude) - 1);
                    int prev_min = 0, prev_max = 0;
                    for (int i = 0; i < sub_dup_nodes.size(); i++) {
                        int start = sub_dup_nodes[i].start;
                        int stop = sub_dup_nodes[i].stop;
                        if (start > prev_min && stop < prev_max) {
                            continue;
                        }
                        if (start < prev_min && stop < prev_max && stop >= prev_min) {
                            stop = prev_min;
                        } else if (start > prev_min && stop > prev_max && start <= prev_max) {
                            start = prev_max;
                        }
                        for (int j = start; j <= stop; j++) {
                            if (duplication_on_path > 0) {
                                many_to_many_labels.insert(j);
                            } else {
                                one_to_many_labels.insert(j);
                            }
                            visited.insert(j);
                        }
                        prev_min = start;
                        prev_max = stop;
                    }

                    for (int i = get<0>(include); i <= get<0>(exclude) - 1; i++) {
                        if (visited.find(i) == visited.end()) {
                            if (duplication_on_path > 0) {
                                one_to_many_labels.insert(i);
                            } else {
                                one_to_one_labels.insert(i);
                            }
                        }
                    }
                }
            } else {
                // if we encounter a dummy duplication label, increment the duplication_on_path counter
                duplication_on_path++;
            }
        }

        for (int label : one_to_one_labels) {
            IndexedGeneTreeNode idx_node = this->gti->leaf_labels.at(label);
            OrthologPair ortho_pair = {
                .gene_name = query_gene,
                .taxon = "",
                .ortholog_name = idx_node.gene_name,
                .type = ONE_TO_ONE
            };
            orthologs.push_back(ortho_pair);
        }
        for (int label : one_to_many_labels) {
            IndexedGeneTreeNode idx_node = this->gti->leaf_labels.at(label);
            OrthologPair ortho_pair = {
                .gene_name = query_gene,
                .taxon = "",
                .ortholog_name = idx_node.gene_name,
                .type = ONE_TO_MANY
            };
            orthologs.push_back(ortho_pair);
        }
        for (int label : many_to_many_labels) {
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

vector<ParalogPair> GeneTree::get_paralogs(string query_gene) {
    vector<ParalogPair> paralogs;
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
        vector<int> paralog_labels;
        for (int j = 0; j < ancestors_idx.size(); j++) {
            int min_label = get<0>(ancestors_idx[j].internal_label);
            int max_label = get<1>(ancestors_idx[j].internal_label);
            if (idx_gene_node.label >= min_label && idx_gene_node.label <= max_label) {
                // label is within the range
                for (int k = min_label; k <= max_label; k++) {
                    if (ancestors_idx[j].node_type == DUPLICATION) {
                            // if the internal node is a speciation node, add all labels in the range that are not visited
                        if (k != idx_gene_node.label && std::find(visited.begin(), visited.end(), k) == visited.end()) {
                            paralog_labels.push_back(k);
                        }
                    }
                    visited.push_back(k);
                }
            }
        }
        wstring tax1 = gene_node->get_taxonomy();
        for (int i = 0; i < paralog_labels.size(); i++) {
            int label = paralog_labels[i];
            IndexedGeneTreeNode idx_node = this->gti->leaf_labels.at(label);
            GeneTreeNode *node = this->root->leaves_map.at(idx_node.node_hash);
            wstring tax2 = node->get_taxonomy();
            if (tax1 == tax2) {
                ParalogPair paralog_pair = {
                    .gene_name = query_gene,
                    .taxon = string(tax1.begin(), tax1.end()),
                    .paralog_name = idx_node.gene_name,
                    .type = WITHIN_SPECIES
                };
                paralogs.push_back(paralog_pair);
            }
        }
    }
    return paralogs;
}

vector<wstring> GeneTree::get_orthologs_naive(string gene_name) {
    vector<wstring> orthologs;
    vector<GeneTreeNode*> genes = this->root->get_leaves();
    for (int i = 0; i < genes.size(); i++) {
        if (genes[i]->get_name() == wstring(gene_name.begin(), gene_name.end())) {
            // found the gene
            GeneTreeNode *gene_node = genes[i];
            vector<GeneTreeNode*> ancestors = gene_node->get_ancestors();
            for (int j = 0; j < ancestors.size(); j++) {
                GeneTreeNode *ancestor = ancestors[j];
                if (ancestor->node_type == DUPLICATION) {
                    // if the ancestor is a duplication node, add all its leaves to the orthologs
                    vector<GeneTreeNode*> leaves = ancestor->get_leaves();
                    for (int k = 0; k < leaves.size(); k++) {
                        if (find(orthologs.begin(), orthologs.end(), leaves[k]->get_name()) == orthologs.end()) {
                            wstring leaf_name = leaves[k]->get_name();
                            orthologs.push_back(leaf_name);
                        }
                    }
                }
            }
        }
    }
    return orthologs;
}