#include "genetree.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include <sstream>
#include <iostream>
#include "wchar.h"
#include "sys/stat.h"

using namespace compara;
using namespace vtdxml;
using namespace std;

GeneTree::GeneTree(const char* filename) {
    this->filename = filename;
    this->parse();
}

GeneTree::~GeneTree() {
    delete this->vn;
    delete this->root;
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

int main (int argc, char *argv[]) {
    GeneTree *gt = new GeneTree(argv[1]);
    gt->print();
    vector<wstring> genes = gt->get_genes();
    for (int i = 0; i < genes.size(); i++) {
        wcout << genes[i] << endl;
    }
    delete gt;
    return 0;
}