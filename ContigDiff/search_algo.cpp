//
// Created by Florian Claisse on 12/05/2023.
//

#include <filesystem>
#include <fstream>
#include <map>
#include <omp.h>
#include "search_algo.h"
#include "../Foundation/include/program_option.h"
#include "../Foundation/include/directory.h"

using namespace std;
namespace fs = std::filesystem;

bool find_inside_file(string &word, unsigned long max_nb_error, ifstream &file) {
    unsigned long nb_error;
    string line;
    while (getline(file, line)) {
        if (line.at(0) == '>') continue;
        for (unsigned long i = 0; i < line.size(); i++) {
            nb_error = 0;
            if ( (i + word.size()) > line.size()) {
                nb_error += ((i + word.size()) - line.size());
                if (nb_error > max_nb_error) break;
            }
            for (unsigned long j = 0; j < word.size(); j++) {
                if (line[i + j] != word[j]) {
                    nb_error++;
                    if (nb_error > max_nb_error) break;
                }
            }
            if (nb_error <= max_nb_error) return true;
        }
    }
    return false;
}

void common_in_two_file(const fs::path &path_A, const fs::path &path_B, const program_option::ContigDiff &options) {
    string filename = directory::fileName(path_A);
    ifstream fileA(path_A);

    vector<map<string, string>> all_common(options.threads);
    string read_in_A, contig_name;
    while(getline(fileA, read_in_A)) {
        // Set contig name
        if (read_in_A.at(0) == '>') { contig_name = read_in_A; continue; }
        // Search
#pragma omp parallel for default(none) shared(read_in_A, path_B, contig_name, options, all_common) num_threads(options.threads)
        for (unsigned long size = read_in_A.size(); size > 0; size--) {
            ifstream input_B(path_B);
            string sub = read_in_A.substr(0, size);
            if (find_inside_file(sub, (size * (100 - options.accept)) / 100, input_B)) {
                all_common[omp_get_thread_num()][contig_name] = sub;
#pragma omp cancel for
            }
            if (input_B.eof()) input_B.clear();
            input_B.seekg(0, ios::beg);
        }
    }
}
