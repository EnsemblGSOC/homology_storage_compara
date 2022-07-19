#include "vtdGen.h"
#include "vtdNav.h"
#include "vtdNav_L5.h"
#include "vtdException.h"
#include "bookMark.h"
#include <iostream>
#include <vector>
#include "genetree_node.h"

namespace compara {
    /**
     * @brief A gene tree object.
     */
    class GeneTree {
        public:
            GeneTree(const char* filename);
            ~GeneTree();
            void parse();
            void print();
            vector<wstring> get_genes();
            vector<wstring> get_orthologs(string gene_name);
            vector<wstring> get_paralogs(string gene_name);
            vtdxml::VTDNav *vn;
        
        private:
            const char* filename;
            compara::GeneTreeNode *root;
            void parse_genetree_node();
    };
}