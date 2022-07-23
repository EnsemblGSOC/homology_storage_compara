#include <iostream>
#include <fstream>
#include <unistd.h>
#include "genetree_index.h"
#include "genetree.h"
#include "time.h"

using namespace std;
using namespace compara;

int main (int argc, char *argv[]) {
    int option;
    char *save_location = NULL;
    char *load_location = NULL;
    char *ortholog_query = NULL;
    char *paralog_query = NULL;
    bool print_tree = false;
    bool print_genes = false;
    bool save_index = false;
    bool load_index = false;
    bool query_orthologs = false;
    bool query_paralogs = false;
    while ((option = getopt(argc, argv, "hpls:i:O:P:")) != -1) {
        switch (option) {
            case 'h':
                cout << "Usage: " << argv[0] << " [OPTIONS] <filename>" << endl;
                cout << "Options:" << endl;
                cout << "  -h\t\t\tPrint this help message" << endl;
                cout << "  -p\t\t\tPrint the gene tree" << endl;
                cout << "  -l\t\t\tList all genes in the gene tree" << endl;
                cout << "  -s\t\t\tConstruct an save index file" << endl;
                cout << "  -i <index_file>\t\tLoad the index file" << endl;
                cout << "  -O <gene name>\t\tList all the orthologs of a gene" << endl;
                cout << "  -P <gene name>\t\tList all the paralogs of a gene" << endl;
                break;
            case 'p':
                print_tree = true;
                break;
            case 'l':
                print_genes = true;
                break;
            case 's':
                save_index = true;
                save_location = optarg;
                break;
            case 'i':
                load_index = true;
                load_location = optarg;
                break;
            case 'O':
                query_orthologs = true;
                ortholog_query = optarg;
                break;
            case 'P':
                query_paralogs = true;
                paralog_query = optarg;
                break;
            default:
                cerr << "Unknown option: " << option << endl;
                return 1;
        }
    }
    if (optind >= argc) {
        cerr << "No input file specified" << endl;
        return 1;
    }
    char* filename = argv[optind];
    
    GeneTree *gt = new GeneTree(filename);
    if (print_tree) {
        gt->print();
    }
    if (print_genes) {
        vector<wstring> genes = gt->get_genes();
        for (auto gene : genes) {
            wcout << gene << endl;
        }
    }
    if (save_index) {
        if (save_location == NULL) {
            cerr << "No index file specified" << endl;
            return 1;
        }
        gt->write_index(save_location);
    }
    if (load_index) {
        if (load_location == NULL) {
            cerr << "No index file specified" << endl;
            return 1;
        }
        gt->load_index(load_location);
    } else {
        gt->write_index("temp.genetreeidx");
        gt->load_index("temp.genetreeidx");
    }
    if (query_orthologs) {
        clock_t c1, c2;
        c1 = clock();
        vector<OrthologPair> genes = gt->get_orthologs(ortholog_query);
        for (auto gene : genes) {
            cout << gene.ortholog_name << " ";
            switch (gene.type) {
                case OrthologType::ONE_TO_ONE:
                    cout << "(1-to-1)";
                    break;
                case OrthologType::ONE_TO_MANY:
                    cout << "(1-to-many)";
                    break;
                case OrthologType::MANY_TO_MANY:
                    cout << "(many-to-many)";
                    break;
                default:
                    cout << "(unknown)";
                    break;
            }
            cout << endl;
        }
        c2 = clock();
        double time_taken = (double)(c2 - c1) / CLOCKS_PER_SEC;
        cout << "Found " << genes.size() << " orthologs in " << time_taken << " seconds" << endl;
    }
    if (query_paralogs) {
        clock_t c1, c2;
        c1 = clock();
        vector<ParalogPair> genes = gt->get_paralogs(paralog_query);
        for (auto gene : genes) {
            cout << gene.paralog_name << " ";
            switch (gene.type) {
                case ParalogType::WITHIN_SPECIES:
                    cout << "(within-species)";
                    break;
                case ParalogType::BETWEEN_SPECIES:
                    cout << "(between-species)";
                    break;
                default:
                    cout << "(unknown)";
                    break;
            }
            cout << endl;
        }
        c2 = clock();
        double time_taken = (double)(c2 - c1) / CLOCKS_PER_SEC;
        cout << "Found " << genes.size() << " paralogs in " << time_taken << " seconds" << endl;
    }
    delete gt;
    return 0;
}