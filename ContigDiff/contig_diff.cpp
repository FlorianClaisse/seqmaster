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
#include "contig_diff.h"

using namespace std;
namespace fs = std::filesystem;

int check_options(const program_option::ContigDiff &options) {
    if (!fs::is_directory(options.inputA)) {
        cout << "Path : " << options.inputA << " n'est pas un dossier.\n";
        return EXIT_FAILURE;
    }

    if (!fs::is_directory(options.inputB)) {
        cout << "Path : " << options.inputB << " n'est pas un dossier.\n";
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}

void find_inside_file(string &word, string &contig_name, unsigned long max_nb_error, ifstream &file, map<string, string> &all_word) {
    string line;
    while(getline(file, line)) {
        if (line.at(0) == '>') continue;
        unsigned long nb_error(0);
        for (unsigned long i = 0; i < line.size(); i++) {
            for (unsigned long j = 0; j < word.size(); j++) {
                if (i + j > line.size()) {
                    nb_error += (word.size() - j);
                    if (nb_error <= max_nb_error) {
                        all_word[word] = contig_name;
                        return;
                    } else break;
                } else if (line[i + j] != word[j]) {
                    nb_error++;
                    if (nb_error > max_nb_error) break;
                }
            }
            if (nb_error <= max_nb_error) {
                all_word[word] = contig_name;
                return;
            }
            nb_error = 0;
        }
    }
}

bool check_inside_file(ifstream &file, const string &word, unsigned long max_nb_error) {
    string line_read;
    unsigned long nb_error(0);
    while(getline(file, line_read)) {
        if (line_read.at(0) == '>') continue;
        for (unsigned long i = 0; i < line_read.size(); i++) {
            for (unsigned long j = 0; j < word.size(); j++) {
                if (i + j > line_read.size()) {
                    nb_error += (word.size() - j);
                    if (nb_error <= max_nb_error) return true;
                    else break;
                } else if (line_read[i + j] != word[j]) {
                    nb_error++;
                    if (nb_error > max_nb_error) break;
                }
            }
            if (nb_error <= max_nb_error) return true;
            nb_error = 0;
        }
    }
    return false;
}

int contig_diff::start(const program_option::ContigDiff &options) {

    if (check_options(options) != EXIT_SUCCESS) return EXIT_FAILURE;

    cout << "Convert input A file's to fastaline.\n";
    unsigned long nb_fasta_files(0);
    for (const auto &currentFile : fs::directory_iterator(options.inputA)) {
        if (!fasta::is_fasta_file(currentFile)) continue;
        if (fasta::is_result_file(currentFile)) continue;
        if (fasta::to_fasta_line(currentFile) != EXIT_SUCCESS) return EXIT_FAILURE;
        nb_fasta_files++;
    }

    cout << "Convert input B file's to fastaline.\n";
    for (const auto &currentFile: fs::directory_iterator(options.inputB)) {
        if (!fasta::is_fasta_file(currentFile)) continue;
        if (fasta::is_result_file(currentFile)) continue;
        if (fasta::to_fasta_line(currentFile) != EXIT_SUCCESS) return EXIT_FAILURE;
    }

    cout << "Check number of files in input A.\n";
    if (nb_fasta_files < 2) {
        cout << "Le nombre de fichier dans l'input A doit être au minimum de 2.\n";
        return EXIT_FAILURE;
    }

    cout << "Find common in 2 A file\n";
    bool first(true);
    ifstream first_input;
    string second_input_path;
    for (const auto &currentFile : fs::directory_iterator(options.inputA)) {
        if (!fasta::is_fastaline_file(currentFile)) continue;
        if (fasta::is_result_file(currentFile)) continue;

        if (first) {
            first = false;
            first_input.open(currentFile.path());
            cout << "\tOpen : " << currentFile.path() << ", in first file.\n";
        } else {
            second_input_path = currentFile.path();
            cout << "\tOpen : " << currentFile.path() << ", in second file.\n";
            break;
        }
    }

    vector<map<string, string>> all_common(options.threads);
    string first_line_read, contig_name;
    while (getline(first_input, first_line_read)) {
        if (first_line_read.at(0) == '>') {
            contig_name = first_line_read;
            continue;
        }
#pragma omp parallel for default(none) shared(first_line_read, all_common, options, second_input_path, contig_name) num_threads(options.threads)
        for (unsigned long size = first_line_read.size(); size > 0; size--) {
            ifstream second_input(second_input_path);
            int currentThread = omp_get_thread_num();
            string sub = first_line_read.substr(0, size);
            if (all_common[currentThread].find(sub) == all_common[currentThread].end()) {
                find_inside_file(sub, contig_name, (size * (100 - options.accept)) / 100, second_input, all_common[currentThread]);
                if (second_input.eof()) second_input.clear();
                second_input.seekg(0, ios::beg);
            }
            second_input.close();
        }
    }

    first_input.close();

    unsigned long common_size(0);
    for (const auto & i : all_common) {
        common_size += i.size();
    }
    cout << "Common start Size : " << common_size << endl;
    cout << "Find common inside all A.\n";

    /*vector<string> key_to_remove;
    for (const auto &currentPath : fs::directory_iterator(options.inputA)) {
        if (!fasta::is_fastaline_file(currentPath)) continue;
        if (fasta::is_result_file(currentPath)) continue;

        ifstream currentFile(currentPath);
        for (const auto &common : all_common) {
            if (!check_inside_file(currentFile, common.first, (common.first.size() * (100 - options.accept)) / 100)) {
                key_to_remove.push_back(common.first);
            }
            if (currentFile.eof()) currentFile.clear();
            currentFile.seekg(0, ios::beg);
        }
        for (const auto &key: key_to_remove) all_common.erase(key);
        key_to_remove.clear();

        currentFile.close();
    }

    cout << "Common end Size : " << all_common.size() << endl;*/


    return EXIT_SUCCESS;
}
