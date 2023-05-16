//
// Created by Florian Claisse on 10/05/2023.
//

#include <iostream>
#include <filesystem>
#include <fstream>
#include <omp.h>
#include <map>
#include <vector>
#include "../Foundation/include/fasta.h"
#include "../Foundation/include/directory.h"
#include "include/contig_diff.h"

// TODO: Optimisation (taille en generale des chaine commune (i++ ou i--))

using namespace std;
namespace fs = std::filesystem;

namespace contig_diff {
    void search_with_2_file(const fs::path &first_path, const fs::path &second_path, map<string, string> &common);
    unsigned long find_inside_file(ifstream &test_file, const string &word);
}

int contig_diff::main(const contig_diff::option &options) {
    cout << "Convert input A file's to fastaline.\n";
    fasta::directory_to_fastaline(options.inputA);

    cout << "Convert input B file's to fastaline.\n";
    fasta::directory_to_fastaline(options.inputB);

    int nb_fasta_file = directory::count_fasta_file(options.inputA);
    if (nb_fasta_file < 2) {
        cout << "Il faut au minimum deux fichiers de type fasta dans le dossier A.\n";
        return EXIT_FAILURE;
    }

    // [contig_value, (contig_name -> filename)]
    map<string, string> common;

    // Init common with 2 file
    pair<fs::path, fs::path> two_first = directory::two_first_fasta(options.inputA);
    search_with_2_file(two_first.first, two_first.second, common);
    cout << common.size() << endl;

    // Compare first file with all files
    for (const auto &path: fs::directory_iterator(options.inputA)) {
        if (path == two_first.first || !directory::is_fastaline_file(path)) continue;
        search_with_2_file(path, two_first.first, common);
        cout << common.size() << endl;
    }

    // Check all possible common inside all file in A.
    unsigned long max_size;
    ifstream *current_file;
    string sub;
    vector<string> value_to_remove;
    map<string, string> value_to_add;
    for (const auto &path: fs::directory_iterator(options.inputA)) {
        if (!directory::is_fastaline_file(path)) continue;
        current_file = directory::read_open(path);
        // Check all common in one file
        for (const auto &value: common) {
            max_size = find_inside_file((*current_file), value.first);
            if (max_size != value.first.size()) {
                sub = value.first.substr(0, max_size);
                if (common.find(sub) == common.end() && value_to_add.find(sub) == value_to_add.end()) {
                    value_to_add[sub] = value.second;
                    value_to_remove.push_back(value.first);
                    cout << "Max size : " << max_size << endl;
                    cout << "Value to remove : " << value.first << endl;
                    cout << "Value to add : " << sub << endl;
                }
            }
        }
        directory::read_close(current_file);
        cout << "Remove : " << value_to_remove.size() << endl;
        cout << "Add : " << value_to_add.size() << endl;

        for (const auto &value: value_to_add) {
            common[value.first] = value.second;
        }
        value_to_add.clear();

        for (const auto &value: value_to_remove) {
            common.erase(value);
        }
        value_to_remove.clear();
    }

    return EXIT_SUCCESS;
}

void contig_diff::search_with_2_file(const fs::path &first_path, const fs::path &second_path, map<string, string> &common) {
    string filename = first_path.stem();
    ifstream *first_file = directory::read_open(first_path);
    ifstream *test_file = directory::read_open(second_path);

    string first_line_read, contig_name;
    while(getline((*first_file), first_line_read)) {
        if (first_line_read.at(0) == '>') {
            contig_name = first_line_read.append(" -> " + filename);
        } else {
            unsigned long size = find_inside_file((*test_file), first_line_read);
            if (size != 0) {
                string sub = first_line_read.substr(0, size);
                 if (common.find(sub) == common.end()) {
                     common[sub] = contig_name;
                 }
            }
            if (test_file->eof()) test_file->clear();
            test_file->seekg(0, ios::beg);
        }
    }

    directory::read_close(test_file);
    directory::read_close(first_file);
}

// Only for 100%
unsigned long contig_diff::find_inside_file(ifstream &test_file, const string &word) {
    unsigned long max_size(0);
    string line_read;
    while(getline(test_file, line_read)) {
        if (line_read.at(0) == '>') continue;
        for (unsigned int i = 1; i <= word.size(); i++) {
            if (line_read.find(word.substr(0, i)) == string::npos) break;
            if (i > max_size) max_size = i;
            if (max_size == word.size()) return word.size();
        }
        if (max_size == word.size()) return word.size();
    }

    return max_size;
}


