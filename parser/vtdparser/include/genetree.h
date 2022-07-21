#include "vtdGen.h"
#include "vtdNav.h"
#include "vtdNav_L5.h"
#include "vtdException.h"
#include "bookMark.h"
#include <iostream>
#include <vector>
#include "genetree_node.h"
#include "genetree_index.h"

namespace compara {
    /**
     * @brief Type of ortholog pair.
     */
    enum {
        ONE_TO_ONE,
        ONE_TO_MANY,
        MANY_TO_MANY
    } OrthologType;

    /**
     * @brief Type of paralog pair.
     */
    enum {
        WITHIN_SPECIES,
        BETWEEN_SPECIES
    } ParalogType;
    
    /**
     * @brief A pair of orthologous genes
     */
    struct OrthologPair {
        string gene_name;
        string taxon;
        string ortholog_name;
        string ortholog_taxon;
    };
    
    /**
     * @brief A gene tree class.
     */
    class GeneTree {
        public:
            GeneTree(const char* filename);
            ~GeneTree();
            void parse();
            void write_index(const char* filename);
            void print();
            void load_index(const char* filename);
            vector<wstring> get_genes();
            vector<string> get_orthologs(string gene_name);
            vector<wstring> get_orthologs_naive(wstring gene_name);
            vector<wstring> get_paralogs(string gene_name);
            vector<wstring> get_paralogs_naive(wstring gene_name);
            vtdxml::VTDNav *vn;
        
        private:
            const char* filename;
            compara::GeneTreeNode *root;
            bool index_loaded;
            GeneTreeIndex *gti;
            void parse_genetree_node();
    };
}