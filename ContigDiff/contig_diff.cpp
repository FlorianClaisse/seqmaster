//
// Created by Florian Claisse on 10/05/2023.
//

#include <iostream>
#include <filesystem>
#include <fstream>
#include <map>
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

bool find_inside_line(string &word, unsigned long pos, string &line, map<string, unsigned long> &all_word) {
    if (all_word.find(word) != all_word.end()) return true;

    if (line.find(word) != string::npos) {
        all_word[word] = pos;
        return true;
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

    cout << "Find common in 2 A file";
    bool first(true);
    ifstream first_input, second_input;
    for (const auto &currentFile : fs::directory_iterator(options.inputA)) {
        if (!fasta::is_fastaline_file(currentFile)) continue;
        if (fasta::is_result_file(currentFile)) continue;

        if (first) {
            first = false;
            first_input.open(currentFile);
            cout << "Open : " << currentFile.path() << ", in first file.\n";
        } else {
            second_input.open(currentFile);
            cout << "Open : " << currentFile.path() << ", in second file.\n";
            break;
        }
    }

    map<string, unsigned long> common;
    string first_line_read, second_line_read, sub;
    // a améliorer
    while (getline(first_input, first_line_read)) {
        if (first_line_read.at(0) == '>') continue;
        cout << "A line read : " << first_line_read << "\n";
        for (unsigned long size = 1; size < first_line_read.size(); size++) {
            for (unsigned long pos = 0; pos <= (first_line_read.size() - size); pos++) {
                sub = first_line_read.substr(pos, size);
                cout << "sub : " << sub << "\n";
                while (getline(second_input, second_line_read)) {
                    if (second_line_read.at(0) == '>') continue;
                    cout << "B line read : " << second_line_read << "\n";
                    if (find_inside_line(sub, pos, second_line_read, common)) {
                        cout << "Find" << "\n";
                        break;
                    }
                }
                if (second_input.eof()) second_input.clear();
                second_input.seekg(0, ios::beg);
            }
        }
    }

    return EXIT_SUCCESS;
}
