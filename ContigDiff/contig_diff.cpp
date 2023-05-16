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

using namespace std;
namespace fs = std::filesystem;

namespace contig_diff {
    void init_with_2_file(const fs::path &first_path, const fs::path &second_path, map<string, string> &common);
    unsigned long find_inside_file(ifstream &test_file, const string &word);
}

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

unsigned long find_inside_file2(ifstream &test_file, const string &word) {
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

void init_with_2_file2(const fs::path &first_path, const fs::path &second_path, map<string, string> &common) {
    string filename = first_path.stem();
    ifstream *first_file = directory::read_open(first_path);

    string first_line_read, contig_name;
    while(getline((*first_file), first_line_read)) {
        if (first_line_read.at(0) == '>') {
            contig_name = first_line_read.append(" -> " + filename);
        } else {
            ifstream *test_file = directory::read_open(second_path);
            unsigned long size = find_inside_file2((*test_file), first_line_read);
            if (size != 0) {
                string sub = first_line_read.substr(0, size);
                if (common.find(sub) == common.end()) {
                    common[sub] = contig_name;
                }
            }
            directory::read_close(test_file);
        }
    }

    directory::read_close(first_file);
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
    map<string, string> common, common2;

    // Init common with 2 file
    pair<fs::path, fs::path> two_first = directory::two_first_fasta(options.inputA);
    init_with_2_file(two_first.first, two_first.second, common);
    //init_with_2_file2(two_first.first, two_first.second, common);
    cout << common.size() << endl;

/*
    vector<map<string, string>> all_contig_common(options.threads);
    string first_line_read, contig_name;
    while (getline(first_input, first_line_read)) {
        if (first_line_read.at(0) == '>') {
            contig_name = first_line_read;
            continue;
        }
#pragma omp parallel for default(none) shared(first_line_read, second_input_path, contig_name, all_contig_common, options) num_threads(options.threads)
        for (unsigned long size = first_line_read.size(); size > 0; size--) {
            ifstream second_input(second_input_path);
            string sub = first_line_read.substr(0, size);
            if (find_inside_file(sub, (size * (100 - options.accept)) / 100, second_input)) {
                all_contig_common[omp_get_thread_num()][contig_name] = sub;
#pragma omp cancel for
            }
            if (second_input.eof()) second_input.clear();
            second_input.seekg(0, ios::beg);
        }
    }

    map<string, string> all_common_result;
    for (const auto &vec : all_contig_common) {
        for (const auto &value: vec) {
            if (all_common_result.find(value.first) == all_common_result.end()) {
                if (value.second > all_common_result[value.first]) all_common_result[value.first] = value.second;
            } else all_common_result[value.first] = value.second;
        }
    }
    all_contig_common.clear();

    cout << "Common start size : " << all_common_result.size() << endl;
*/
/*
    all_common.clear();

    cout << "Find common inside all A.\n";
    long nb_files = distance(fs::directory_iterator(options.inputA), fs::directory_iterator{}); // C++17

    vector<vector<string>> all_key_to_remove(options.threads);
#pragma omp parallel for default(none) shared(nb_files, options, common_result, all_key_to_remove) num_threads(options.threads)
    for (long i = 0; i < nb_files; i++) {
        fs::directory_entry currentPath = *next(fs::directory_iterator(options.inputA), i);
        if (!fasta::is_fastaline_file(currentPath)) continue;
        if (fasta::is_result_file(currentPath)) continue;

        ifstream currentFile(currentPath);
        for (const auto &common : common_result) {
            if (!check_inside_file(currentFile, common.first, (common.first.size() * (100 - options.accept)) / 100)) {
                all_key_to_remove[omp_get_thread_num()].push_back(common.first);
            }
            if (currentFile.eof()) currentFile.clear();
            currentFile.seekg(0, ios::beg);
        }

        currentFile.close();
    }

    // TODO: Garder que le plus grande chaine commune pour chaque contig

    // Remove all keys
    for (const auto &vec: all_key_to_remove) {
        for (const auto &key: vec) {
            common_result.erase(key);
        }
    }
    all_key_to_remove.clear();

    cout << "Common end Size : " << common_result.size() << endl;

    cout << "Find inside B.\n";
    nb_files = distance(fs::directory_iterator(options.inputB), fs::directory_iterator{}); // C++17

    vector<vector<string>> all_key_to_remove_in_B(options.threads);
#pragma omp parallel for default(none) shared(nb_files, options, common_result, all_key_to_remove_in_B) num_threads(options.threads)
    for (long i = 0; i < nb_files; i++) {
        fs::directory_entry currentPath = *next(fs::directory_iterator(options.inputB), i);
        if (!fasta::is_fastaline_file(currentPath)) continue;
        if (fasta::is_result_file(currentPath)) continue;

        ifstream currentFile(currentPath);
        for (const auto &common : common_result) {
            if (check_inside_file(currentFile, common.first, (common.first.size() * (100 - options.accept)) / 100)) {
                all_key_to_remove_in_B[omp_get_thread_num()].push_back(common.first);
            }
            if (currentFile.eof()) currentFile.clear();
            currentFile.seekg(0, ios::beg);
        }

        currentFile.close();
    }

    // Remove all keys
    for (const auto &vec: all_key_to_remove_in_B) {
        for (const auto &key: vec) {
            common_result.erase(key);
        }
    }
    all_key_to_remove_in_B.clear();

    cout << "Common end Size : " << common_result.size() << endl;*/

    return EXIT_SUCCESS;
}

void contig_diff::init_with_2_file(const fs::path &first_path, const fs::path &second_path, map<string, string> &common) {
    string filename = first_path.stem();
    ifstream *first_file = directory::read_open(first_path);

    string first_line_read, contig_name;
    while(getline((*first_file), first_line_read)) {
        if (first_line_read.at(0) == '>') {
            contig_name = first_line_read.append(" -> " + filename);
        } else {
            ifstream *test_file = directory::read_open(second_path);
            unsigned long size = find_inside_file((*test_file), first_line_read);
            if (size != 0) {
                string sub = first_line_read.substr(0, size);
                 if (common.find(sub) == common.end()) {
                     common[sub] = contig_name;
                 }
            }
            directory::read_close(test_file);
        }
    }

    directory::read_close(first_file);
}

// Only for 100%
unsigned long contig_diff::find_inside_file(ifstream &test_file, const string &word) {
    unsigned long max_size(0);
    bool error_found(false);
    string line_read;
    while(getline(test_file, line_read)) {
        if (line_read.at(0) == '>') continue;
        for (unsigned long offset = 0; offset < line_read.size(); offset++) {
            for (unsigned long i = 0; i < word.size(); i++) {
                word.find("test");
                if (offset + i > line_read.size()) {
                    if (i > max_size) max_size = i;
                    error_found = true;
                    break;
                } else if (line_read[offset + i] != word[i]) {
                    if ((i > max_size)) max_size = i;
                    error_found = true;
                    break;
                }
            }
            if (!error_found) return word.size();
            error_found = false;
        }
    }
    return max_size;
}


