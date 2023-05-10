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

bool find_inside_file(string &word, string &contig_name, unsigned long max_nb_error, ifstream &file, map<string, string> &all_word) {
    if (all_word.find(word) != all_word.end()) return true;

    unsigned long nb_error(0);
    string line;
    while(getline(file, line)) {
        if (line.at(0) == '>') continue;
        for (unsigned long i = 0; i < line.size(); i++) {
            for (unsigned long j = 0; j < word.size(); j++) {
                if (i + j > line.size()) {
                    nb_error += (word.size() - j);
                    break;
                }
                if (line[i + j] != word[j]) nb_error++;
                if (nb_error > max_nb_error) break;
            }
            if (nb_error <= max_nb_error) {
                all_word[word] = contig_name;
                if (file.eof()) file.clear();
                file.seekg(0, ios::beg);
                return true;
            }
            nb_error = 0;
        }
    }

    if (file.eof()) file.clear();
    file.seekg(0, ios::beg);
    return false;
}

// TODO: Continuer cette fonction.
void check_inside_file(ifstream &file, map<string, string> &all_word) {
    unsigned long nb_error(0);
    string line_read;
    while(getline(file, line_read)) {
        for (const auto &word : all_word) {
            for (unsigned long i = 0; i < line_read.size(); i++) {
                for (unsigned long j = 0; j < word.first.size(); j++) {
                    if (i + j > line_read.size()) {
                        nb_error += (word.size() - j);
                        break;
                    }
                    if (line[i + j] != word[j]) nb_error++;
                    if (nb_error > max_nb_error) break;
                }
                if (nb_error <= max_nb_error) {
                    all_word[word] = contig_name;
                    if (file.eof()) file.clear();
                    file.seekg(0, ios::beg);
                    return true;
                }
                nb_error = 0;
            }
        }
    }
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

    map<string, string> all_common;
    string first_line_read, second_line_read, sub, contig_name;
    while (getline(first_input, first_line_read)) {
        if (first_line_read.at(0) == '>') {
            contig_name = first_line_read;
            continue;
        }
        for (unsigned long size = first_line_read.size(); size > 0; size--) {
            sub = first_line_read.substr(0, size);
            find_inside_file(sub, contig_name, (size * (100 - options.accept)) / 100, second_input, all_common);
        }
    }

    cout << "Common Size : " << all_common.size() << endl;

    for (const auto &currentPath : fs::directory_iterator(options.inputA)) {
        if (!fasta::is_fastaline_file(currentPath)) continue;
        if (fasta::is_result_file(currentPath)) continue;

        ifstream currentFile(currentPath);
        check_inside_file(currentFile, all_common);
    }


    return EXIT_SUCCESS;
}
