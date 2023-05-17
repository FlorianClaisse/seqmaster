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

// TODO: Optimisation
// TODO: Amelioration supprimer les chaines doubles

using namespace std;
namespace fs = std::filesystem;

namespace contig_diff {
    typedef struct {
        string filename;
        string inputB;
        int percentage;
    } Common;

    void search_with_2_file(const fs::path &first_path, const fs::path &second_path, map <string, Common> &common);
    unsigned long find_inside_file(ifstream &test_file, const string &word);
    void check_on_directory(const fs::path &directory_path, map <string, Common> &common);
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

    // [sub_contig_value, Common]
    map<string, Common> common;

    // Init common with 2 file
    pair<fs::path, fs::path> two_first = directory::two_first_fasta(options.inputA);
    search_with_2_file(two_first.first, two_first.second, common);

    // Compare first file with all files
    for (const auto &path: fs::directory_iterator(options.inputA)) {
        if (path == two_first.first || !directory::is_fastaline_file(path)) continue;
        search_with_2_file(path, two_first.first, common);
        cout << common.size() << endl;
    }

    // Check all possible common inside all file in A.
    check_on_directory(options.inputA, common);


    return EXIT_SUCCESS;
}

void contig_diff::search_with_2_file(const fs::path &first_path, const fs::path &second_path, map<string, Common> &common) {
    string filename = first_path.stem();
    ifstream *first_file = directory::read_open(first_path);
    ifstream *test_file = directory::read_open(second_path);

    string first_line_read, contig_name;
    while(getline((*first_file), first_line_read)) {
        if (first_line_read.at(0) == '>') {
            contig_name = first_line_read.append(" -> " + filename);
            contig_name.erase(remove(contig_name.begin(), contig_name.end(), '\n'), contig_name.end());
            contig_name.erase(remove(contig_name.begin(), contig_name.end(), '\r'), contig_name.end());
        } else {
            unsigned long size = find_inside_file((*test_file), first_line_read);
            if (size != 0) {
                string sub = first_line_read.substr(0, size);
                 if (common.find(sub) == common.end()) {
                     common[sub] = {contig_name, "", 100};
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

void contig_diff::check_on_directory(const fs::path &directory_path, map<string, Common> &common) {
    unsigned long max_size;
    ifstream *current_file;
    string sub;
    vector<string> value_to_remove;
    map<string, Common> value_to_add;
    for (const auto &path: fs::directory_iterator(directory_path)) {
        if (!directory::is_fastaline_file(path)) continue;
        current_file = directory::read_open(path);
        // Check all common in one file
        for (const auto &value: common) {
            max_size = find_inside_file((*current_file), value.first);
            if (max_size != value.first.size()) {
                sub = value.first.substr(0, max_size);
                if ( max_size != 0) {
                    value_to_add[value.first] = {sub, "", 100};
                }
                value_to_remove.push_back(value.first);
                cout << "Max size : " << max_size << ", Value to remove : " << value.first << ", Value to add : " << sub << endl;
            }

            if (current_file->eof()) current_file->clear();
            current_file->seekg(0, ios::beg);
        }
        directory::read_close(current_file);
        cout << "Remove : " << value_to_remove.size() << ", Add : " << value_to_add.size() << endl;

        for (const auto &value: value_to_remove) {
            common.erase(value);
        }
        value_to_remove.clear();

        for (const auto &value: value_to_add) {
            common[value.first] = value.second;
        }
        value_to_add.clear();
    }
}



