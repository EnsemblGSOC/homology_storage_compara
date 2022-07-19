#include "genetree.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include <sstream>
#include <iostream>
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
}

GeneTreeNode::~GeneTreeNode() {
    delete this->bm;
    for (vector<GeneTreeNode*>::iterator it = this->children.begin(); it != this->children.end(); ++it) {
        delete *it;
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
        // inherit the ancestors of the parent and add this node to the ancestors
        this->children[i]->ancestors = this->ancestors;
        this->children[i]->ancestors.push_back(this);
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
    }
    if (this->is_leaf())
        this->node_type = LEAF;
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

/**
 * @brief Get the ancestors of the current node.
 * 
 * @return vector<GeneTreeNode*> 
 */
vector<GeneTreeNode*> GeneTreeNode::get_ancestors() {
    return this->ancestors;
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
            wcout << "(specation) ";
        }
        if (children[i]->node_type == DUPLICATION) {
            wcout << "(duplication) ";
        }
        wcout << children[i]->get_name() << endl;
        children[i]->print(depth + 1);
    }
}