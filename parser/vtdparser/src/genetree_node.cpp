#include "indexed_genetree_node.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include <sstream>
#include <iostream>
#include <fstream>
#include <clocale>
#include "wchar.h"
#include "sys/stat.h"

using namespace compara;
using namespace vtdxml;
using namespace std;

static inline std::wstring string_lower(const std::wstring &str) {
    std::wstring strcopy(str.size(), 0);
    std::transform(str.begin(), str.end(), strcopy.begin(), towlower);
    return strcopy;
}

GeneTreeNode::GeneTreeNode(BookMark *bm) {
    this->bm = bm;
    this->load_children();
    this->load_node_type();
    this->parent = NULL;
}

GeneTreeNode::~GeneTreeNode() {
    delete this->bm;
    for (int i = 0; i < this->children.size(); i++) {
        if (children[i] != NULL) {
            delete children[i];
        }
    }
}

/**
 * @brief Get the first occurrence of a node.
 * 
 * @param node_name 
 * @return wstring 
 */
wstring GeneTreeNode::get_first_node(wstring node_name, wstring attrib_name) {
    VTDNav *vn = this->bm->getNav();
    // go to the next level
    vn->toElement(FIRST_CHILD, const_cast<wchar_t *>(node_name.c_str()));
    int name_index = vn->getCurrentIndex();
    if (name_index != -1) {
        wstring name_str = vn->toStringLowerCase(name_index);
        if (name_index != -1) {
            if (vn->matchElement(node_name.c_str())) {
                if (attrib_name.empty()) {
                    if (vn->getText() != -1)
                        return vn->toString(vn->getText());
                } else {
                    if (vn->getAttrVal(attrib_name.c_str()) != -1)
                        return vn->toString(vn->getAttrVal(attrib_name.c_str()));
                }
            }
        }
    }
    // return the cursor to the parent
    vn->toElement(PARENT);
    return L"";
}

/**
 * @brief Get the gene name of the node.
 * 
 * @return wstring - gene name, empty string if the node is not a gene
 */
wstring GeneTreeNode::get_name() {
    return this->get_first_node(L"name");
}

/**
 * @brief Get the type of the current node.
 * 
 * @return wstring - type of node
 */
wstring GeneTreeNode::get_element_type() {
    VTDNav *vn = this->bm->getNav();
    int idx = vn->getCurrentIndex();
    if (idx != -1)
        return vn->toStringLowerCase(idx);
    else
        return L"";
}

/**
 * @brief Get the children of this node.
 * 
 * @return vector<GeneTreeNode*> 
 */
vector<GeneTreeNode*> GeneTreeNode:: get_children() {
    return this->children;
}

/**
 * @brief Load the children of the subtree rooted at this node recursively.
 * 'leaves_map' is also populated as we backtrace from the recursion.
 * 
 */
void GeneTreeNode::load_children() {
    VTDNav *vn = this->bm->getNav();
    vector<GeneTreeNode*> children;
    // go to the next level
    vn->toElement(FIRST_CHILD);
    int child_index = vn->getCurrentIndex();
    if (child_index != -1) {
        wstring child_str = L"";
        // go through the siblings until we are at a clade node
        while (child_index != -1) {
            child_index = vn->getCurrentIndex();
            child_str = vn->toStringLowerCase(child_index);
            // if we find a clade node, create a new GeneTreeNode and add it to the vector
            if (child_str.find(L"clade") != wstring::npos) {
                VTDNav *vn_cpy = vn->cloneNav();
                // TODO: use smart pointer
                GeneTreeNode *child = new GeneTreeNode(new BookMark(vn_cpy));
                children.push_back(child);
                if (child->is_leaf()) {
                    this->leaves_map.insert(pair<int, GeneTreeNode*>(child->bm->getNav()->hashCode(), child));
                }
            }
            int jump_result = vn->toElement(NEXT_SIBLING);
            if (jump_result == 0) {
                break;
            }
        }
    }
    // return the cursor to the parent
    vn->toElement(PARENT);
    this->children = children;
    for (int i = 0; i < this->children.size(); i++) {
        this->children[i]->parent = this;
        // this does not preserve entries in the source
        // in the end, only the root will have a non-empty leaves_map
        this->leaves_map.merge(this->children[i]->leaves_map);
    }
}

void GeneTreeNode::load_node_type() {
    VTDNav *vn = this->bm->getNav();
    if (vn->toElement(FIRST_CHILD, const_cast<wchar_t*>(L"events"))) {
        if (vn->toElement(FIRST_CHILD, const_cast<wchar_t*>(L"speciations"))) {
            if (vn->getCurrentIndex() != -1) {
                this->node_type = SPECIATION;
            }
            vn->toElement(PARENT);
        }
        if (vn->toElement(FIRST_CHILD, const_cast<wchar_t*>(L"duplications"))) {
            if (vn->getCurrentIndex() != -1) {
                this->node_type = DUPLICATION;
            }
            vn->toElement(PARENT);
        }
        vn->toElement(PARENT);
    } else {
        if (this->is_leaf())
            this->node_type = LEAF;
        else
            this->node_type = OTHER;
    }
    if (this->node_type == DUPLICATION) {
        if (this->get_confidence_score() <= 0) {
            this->node_type = DUBIOUS;
        }
    }
}

double GeneTreeNode::get_confidence_score() {
    VTDNav *vn = this->bm->getNav();
    if (vn->toElement(LAST_CHILD, const_cast<wchar_t*>(L"confidence"))) {
        if (vn->hasAttr(L"type")) {
            wstring type = vn->toString(vn->getAttrVal(L"type"));
            if (type.compare(L"duplication_confidence_score") == 0) {
                if (vn->getCurrentIndex() != -1) {
                    wstring confidence_str = vn->toString(vn->getText());
                    return stod(confidence_str);
                }
            }
        }
        vn->toElement(PARENT);
    }
    return 0.0;
}

/**
 * @brief Check if the node is a leaf.
 * 
 * @return true 
 * @return false 
 */
bool GeneTreeNode::is_leaf() {
    return this->children.size() == 0;
}

/**
 * @brief Check if the node is a root.
 * 
 * @return true 
 * @return false 
 */
bool GeneTreeNode::is_root() {
    return this->parent == NULL;
}

/**
 * @brief Get the descendants of the current node using a depth-first search.
 * 
 * @return vector<GeneTreeNode*> 
 */
vector<GeneTreeNode*> GeneTreeNode::get_descendants() {
    vector<GeneTreeNode*> descendants;
    for (int i = 0; i < this->children.size(); i++) {
        vector<GeneTreeNode*> child_descendants = this->children[i]->get_descendants();
        for (int j = 0; j < child_descendants.size(); j++) {
            descendants.push_back(child_descendants[j]);
        }
        descendants.push_back(this->children[i]);
    }
    return descendants;
}

vector<GeneTreeNode*> GeneTreeNode::get_leaves() {
    vector<GeneTreeNode*> leaves;
    vector<GeneTreeNode*> descendants = this->get_descendants();
    for (int i = 0; i < descendants.size(); i++) {
        if (descendants[i]->is_leaf()) {
            leaves.push_back(descendants[i]);
        }
    }
    return leaves;
}

/**
 * @brief Get the ancestors of the current node.
 * 
 * @return vector<GeneTreeNode*> 
 */
vector<GeneTreeNode*> GeneTreeNode::get_ancestors() {
    vector <GeneTreeNode*> ancestors;
    GeneTreeNode *curr = this;
    while (curr != NULL && curr->parent != NULL) {
        ancestors.push_back(curr->parent);
        curr = curr->parent;
    }
    return ancestors;
}

/**
 * @brief Print the gene tree rooted at this node.
 * 
 * @param depth indentation
 */
void GeneTreeNode::print(int depth) {
    vector<GeneTreeNode*> children;
    children = this->get_children();
    for (int i = 0; i < children.size(); i++) {
        for (int j = 0; j < depth; j++) {
            cout << "\t";
        }
        wcout << children[i]->get_element_type();
        wcout << " ";
        if (children[i]->node_type == SPECIATION) {
            wcout << "(speciation) ";
        }
        if (children[i]->node_type == DUPLICATION) {
            wcout << "(duplication" << " " << children[i]->get_confidence_score() << ") ";
        }
        wcout << children[i]->get_name() << endl;
        children[i]->print(depth + 1);
    }
}

void GeneTreeNode::write_index(int label, ostream &out) {
    VTDNav *vn = this->bm->getNav();
    if (this->node_type == LEAF) {
        wstring name = this->get_name();
        string narrowed_name = string(name.begin(), name.end());
        IndexedGeneTreeNode indexed_node = IndexedGeneTreeNode(this->bm->hashCode(), label, narrowed_name);
        indexed_node.write_index(out);
    }
}

void GeneTreeNode::write_index(tuple<int, int> label, ostream &out) {
    IndexedGeneTreeNode indexed_node = IndexedGeneTreeNode(this->bm->hashCode(), label, this->node_type);
    indexed_node.write_index(out);
}