#include <iostream>
#include <fstream>
#include <unistd.h>
#include "genetree_index.h"
#include "genetree.h"
#include "time.h"
#include "matplotlibcpp.h"
#include <filesystem>

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

int main (int argc, char *argv[]) {
    clock_t c1, c2;
    char *dir = argv[1];

    vector<int> sizes;
    vector<double> speeds;
    vector<double> times;
    vector<int> ortho_counts;
    vector<int> heights;
    vector<int> total_nodes;
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
    plt::save("cpp_benchmark.pdf");
    plt::show();
}