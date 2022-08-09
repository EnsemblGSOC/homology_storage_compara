#include <iostream>
#include <fstream>
#include <unistd.h>
#include "genetree_index.h"
#include "genetree.h"
#include "time.h"
#include "matplotlibcpp.h"
#include <filesystem>
#include <libgen.h>

using namespace std;
using namespace compara;
namespace plt = matplotlibcpp;

int randint (int n) {
    if ((n - 1) == RAND_MAX) {
        return rand();
    } else {
        // Supporting larger values for n would requires an even more
        // elaborate implementation that combines multiple calls to rand()
        assert (n <= RAND_MAX);
        // Chop off all of the values that would cause skew...
        int end = RAND_MAX / n; // truncate skew
        assert (end > 0);
        end *= n;
        // ... and ignore results from rand() that fall above that limit.
        // (Worst case the loop condition should succeed 50% of the time,
        // so we can expect to bail out of this loop pretty quickly.)
        int r;
        while ((r = rand()) >= end);
        return r % n;
    }
}

void benchmark_query(char *dir) {
    clock_t c1, c2;

    vector<int> sizes;
    vector<double> speeds;
    vector<double> times;
    vector<int> ortho_counts;
    vector<int> heights;
    vector<int> total_nodes;
    vector<int> dup_counts;
    vector<int> xs;

    for (const auto & entry : filesystem::directory_iterator(dir)) {
        string filename_str = entry.path().string();
        const char *filename = filename_str.c_str();
        GeneTree *gt = new GeneTree(filename);
        vector<wstring> leaves = gt->get_genes();
        int size = leaves.size();
        // if (size < 100) {
        //     continue;
        // }
        gt->write_index("temp.genetreeidx");
        gt->load_index("temp.genetreeidx");
        c1 = clock();
        int ortho_count = 0;
        for (int i = 0; i < 100000; i++) {
            int rand_index = randint(size);
            // int rand_index = size / 2;
            wstring name = leaves[rand_index];
            ortho_count += gt->get_orthologs(string(name.begin(), name.end())).size();
        }
        c2 = clock();

        double speed = 100000 / ((double)(c2 - c1) / CLOCKS_PER_SEC);

        int height = gt->root->get_height();
        vector<GeneTreeNode*> nodes = gt->root->get_descendants();
        int total = gt->root->get_descendants().size();
        int total_dup = 0;

        for (int i = 0; i < nodes.size(); i++) {
            if (nodes[i]->node_type == DUPLICATION)
                total_dup++;
        }

        sizes.push_back(size);
        heights.push_back(height);
        ortho_counts.push_back(ortho_count / 100000);
        total_nodes.push_back(total);
        dup_counts.push_back(total_dup);

        xs.push_back(height * log10(total_dup) + ortho_count / 100000);
        
        speeds.push_back(100000 / ((c2 - c1) / (double)CLOCKS_PER_SEC));
        times.push_back((c2 - c1) / (double)CLOCKS_PER_SEC / 100000);

        cout << "Tree size: " << size << endl;
        cout << "File: " << filename << endl;
        cout << "Avg ortholog count: " << ortho_count / 100000 << endl;
        cout << "Time for 100000 ortholog queries: " << (c2 - c1) / (double)CLOCKS_PER_SEC << endl;
        cout << "Speed: " << 100000 / ((c2 - c1) / (double)CLOCKS_PER_SEC) << " queries/sec" << endl;
        cout << "----------------------------------------------------" << endl;
    }
    plt::plot(xs, times, "ob");
    plt::xlabel("h * log(d) + o");
    plt::ylabel("Time per ortholog query (s)");
    plt::save("cpp_benchmark_query.pdf");
    plt::show();
    ofstream outfile("cpp_benchmark_query.csv");
    outfile << "size,speed,time,ortholog_count,height,total_nodes,dup_nodes" << endl;
    for (int i = 0; i < sizes.size(); i++) {
        outfile << sizes[i] << "," << speeds[i] << "," << times[i] << "," << ortho_counts[i] << "," << heights[i] << "," << total_nodes[i] << "," << dup_counts[i] << endl;
    }
    outfile.close();
}

void benchmark_loading(char *dir, char *idx_dir) {
    clock_t c1, c2;
    vector<double> speeds;
    vector<double> times;
    vector<int> sizes;
    vector<int> total_nodes;
    vector<int> xs;
    vector<double> file_sizes;
    vector<int> dup_counts;

    for (const auto & entry : filesystem::directory_iterator(dir)) {
        string filename = entry.path().string();
        string base_filename = filename.substr(filename.find_last_of("/\\") + 1);
        const char *filename_c = filename.c_str();
        // find corresponding index file
        string idx_filename = string(idx_dir) + "/" + base_filename + ".gtidx";
        const char *idx_filename_c = idx_filename.c_str();
        c1 = clock();
        for (int i = 0; i < 100; i++) {
            GeneTree *gt = new GeneTree(filename_c);
            gt->load_index(idx_filename_c);
        }
        c2 = clock();
        double speed = 100 / ((double)(c2 - c1) / CLOCKS_PER_SEC);
        double time = (c2 - c1) / (double)CLOCKS_PER_SEC;
        speeds.push_back(speed);
        times.push_back(time / 100);

        GeneTree *gt = new GeneTree(filename_c);
        gt->load_index(idx_filename_c);
        gt = new GeneTree(filename_c);
        gt->load_index(idx_filename_c);
        int size = gt->get_genes().size();
        sizes.push_back(size);

        vector<GeneTreeNode*> nodes = gt->root->get_descendants();
        int total = nodes.size();
        int total_dup = 0;

        for (int i = 0; i < total; i++) {
            if (nodes[i]->node_type == DUPLICATION)
                total_dup++;
        }
        int x = size + total + total_dup;
        xs.push_back(x);

        dup_counts.push_back(total_dup);
        total_nodes.push_back(total);

        double filesize = filesystem::file_size(entry.path()) + filesystem::file_size(idx_filename);
        // filesize in MB
        filesize /= (1024 * 1024);
        file_sizes.push_back(filesize);

        cout << "File: " << filename << " with " << size << " genes" << endl;
        cout << "Index file: " << idx_filename << endl;
        cout << "Total size of tree and index: " << filesize << " MB" << endl;
        cout << "Time for 100 loading: " << time << endl;
        cout << "Speed: " << speed << " trees/sec" << endl;
        cout << "Speed: " << speed * size << " genes/sec" << endl;
        cout << "Speed: " << filesize / (time / 100) << " MB/s" << endl;
        cout << "----------------------------------------------------" << endl;
    }
    plt::figure(1);
    plt::plot(file_sizes, times, "ob");
    plt::xlabel("file size (MB)");
    plt::ylabel("Time per loading (s)");
    plt::save("cpp_benchmark_loading.pdf");
    plt::figure(2);
    plt::plot(sizes, times, "ob");
    plt::xlabel("number of genes");
    plt::ylabel("Time per loading (s)");
    plt::save("cpp_benchmark_loading2.pdf");
    
    ofstream outfile("cpp_benchmark_loading.csv");
    outfile << "size,speed,time,file_size,total_nodes,dup_nodes" << endl;
    for (int i = 0; i < sizes.size(); i++) {
        outfile << sizes[i] << "," << speeds[i] << "," << times[i] << "," << file_sizes[i] << "," << total_nodes[i] << "," << dup_counts[i] << endl;
    }
    outfile.close();
}

void benchmark_combined(char *dir) {
    clock_t c1, c2;

    vector<int> sizes;
    vector<double> speeds;
    vector<double> times;
    vector<int> ortho_counts;
    vector<int> heights;
    vector<int> total_nodes;
    vector<int> dup_counts;
    vector<int> xs;

    for (const auto & entry : filesystem::directory_iterator(dir)) {
        string filename_str = entry.path().string();
        const char *filename = filename_str.c_str();

        GeneTree *gt = new GeneTree(filename);
        vector<wstring> leaves = gt->get_genes();
        int size = leaves.size();

        gt->write_index("temp.genetreeidx");

        c1 = clock();
        int ortho_count = 0;
        for (int i = 0; i < 1000; i++) {
            GeneTree *gt = new GeneTree(filename);
            gt->load_index("temp.genetreeidx");
            int rand_index = randint(size);
            // int rand_index = size / 2;
            wstring name = leaves[rand_index];
            ortho_count += gt->get_orthologs(string(name.begin(), name.end())).size();
        }
        c2 = clock();

        double speed = 1000 / ((double)(c2 - c1) / CLOCKS_PER_SEC);

        int height = gt->root->get_height();
        vector<GeneTreeNode*> nodes = gt->root->get_descendants();
        int total = gt->root->get_descendants().size();
        int total_dup = 0;

        for (int i = 0; i < nodes.size(); i++) {
            if (nodes[i]->node_type == DUPLICATION)
                total_dup++;
        }

        sizes.push_back(size);
        heights.push_back(height);
        ortho_counts.push_back(ortho_count / 1000);
        total_nodes.push_back(total);
        dup_counts.push_back(total_dup);

        xs.push_back(height * log10(total_dup) + ortho_count / 1000);
        
        speeds.push_back(1000 / ((c2 - c1) / (double)CLOCKS_PER_SEC));
        times.push_back((c2 - c1) / (double)CLOCKS_PER_SEC / 1000);

        cout << "Tree size: " << size << endl;
        cout << "File: " << filename << endl;
        cout << "Avg ortholog count: " << ortho_count / 1000 << endl;
        cout << "Time for 1000 ortholog queries: " << (c2 - c1) / (double)CLOCKS_PER_SEC << endl;
        cout << "Time per ortholog query: " << (c2 - c1) / (double)CLOCKS_PER_SEC / 1000 << endl;
        cout << "Speed: " << 1000 / ((c2 - c1) / (double)CLOCKS_PER_SEC) << " queries/sec" << endl;
        cout << "----------------------------------------------------" << endl;
    }
    plt::plot(sizes, times, "ob");
    plt::xlabel("number of genes");
    plt::ylabel("Time per ortholog query (s)");
    plt::save("cpp_benchmark_combined.pdf");
    plt::show();
    ofstream outfile("cpp_benchmark_combined.csv");
    outfile << "size,speed,time,ortholog_count,height,total_nodes,dup_nodes" << endl;
    for (int i = 0; i < sizes.size(); i++) {
        outfile << sizes[i] << "," << speeds[i] << "," << times[i] << "," << ortho_counts[i] << "," << heights[i] << "," << total_nodes[i] << "," << dup_counts[i] << endl;
    }
    outfile.close();
}

int main (int argc, char *argv[]) {
    int option;
    while ((option = getopt(argc, argv, "qol:")) != -1) {
        switch (option) {
            case 'q':
                benchmark_query(argv[optind]);
                break;
            case 'l':
                benchmark_loading(argv[optind], optarg);
                break;
            case 'o':
                benchmark_combined(argv[optind]);
                break;
            default:
                cout << "Usage: " << argv[0] << " [-q] [-l <index dir>] [-o]" << endl;
                cout << " -q: benchmark query" << endl;
                cout << " -l: benchmark loading" << endl;
                cout << " -o: benchmark combined (loading + query)" << endl;
                break;
        }
    }
}