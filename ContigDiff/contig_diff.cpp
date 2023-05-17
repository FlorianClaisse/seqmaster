//
// Created by Florian Claisse on 10/05/2023.
//

#include <iostream>
#include <filesystem>
#include <fstream>
#include <map>
#include <vector>
#include <algorithm>
#include "../Utils/include/fasta.h"
#include "../Utils/include/directory.h"
#include "include/contig_diff.h"

// TODO: Optimisation: MultiThreads
// TODO: Amelioreation: Prendre en compte un pourcentage different de 100%

using namespace std;
namespace fs = std::filesystem;

namespace contig_diff {
    typedef struct {
        string filename;
        string contigName;
        string inputB;
        int percentage;
    } Common;

    void search_with_2_file(const fs::path &first_path, const fs::path &second_path, map <string, Common> &common);
    unsigned long find_inside_file(ifstream &test_file, const string &word);
    void check_on_directory(const fs::path &directory_path, map<string, Common> &common);
    void verif_on_directory(const fs::path &directory_path, map<string, Common> &common);
    void generate_output(const fs::path &directory_path, const map<string, Common> &common);
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

    cout << "Start finding all common inside A.\n";

    // Init common with 2 file
    pair<fs::path, fs::path> two_first = directory::two_first_fasta(options.inputA);
    search_with_2_file(two_first.first, two_first.second, common);

    // Compare first file with all files
    for (const auto &path: fs::directory_iterator(options.inputA)) {
        if (path == two_first.first || !directory::is_fastaline_file(path)) continue;
        search_with_2_file(path, two_first.first, common);
    }

    // Check all possible common inside all file in A.
    check_on_directory(options.inputA, common);
    cout << "Finish\nAll common inside A size : " << common.size() << endl;
    cout << "Start check all common in A with B.\n";
    // Check all common in A with B;
    verif_on_directory(options.inputB, common);
    cout << "Finish\nAll common after check in B size : " << common.size() << endl;
    cout << "Start Generate output\n";
    generate_output(options.output, common);
    cout << "Finish\n";

    return EXIT_SUCCESS;
}

void contig_diff::search_with_2_file(const fs::path &first_path, const fs::path &second_path, map<string, Common> &common) {
    string filename = first_path.stem();
    ifstream *first_file = directory::read_open(first_path);
    ifstream *test_file = directory::read_open(second_path);

    string first_line_read, contig_name;
    while(getline((*first_file), first_line_read)) {
        if (first_line_read.at(0) == '>') {
            contig_name = first_line_read;
            contig_name.erase(remove_if(contig_name.begin(), contig_name.end(), [](char c) { return c == '\n' || c == '\r'; }), contig_name.end());
        } else {
            unsigned long size = find_inside_file((*test_file), first_line_read);
            if (size != 0) {
                string sub = first_line_read.substr(0, size);
                 if (common.find(sub) == common.end()) {
                     common[sub] = {filename, contig_name, "", 100};
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
                if (common.find(sub) == common.end() && value_to_add.find(sub) == value_to_add.end() && max_size != 0) {
                    value_to_add[sub] = {value.second.filename, value.second.contigName, "", 100};
                }
                value_to_remove.push_back(value.first);
            }

            if (current_file->eof()) current_file->clear();
            current_file->seekg(0, ios::beg);
        }
        directory::read_close(current_file);

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

void contig_diff::verif_on_directory(const fs::path &directory_path, map<std::string, contig_diff::Common> &common) {
    unsigned long max_size;
    ifstream *current_file;
    string sub;
    for (const auto &path: fs::directory_iterator(directory_path)) {
        if (!directory::is_fastaline_file(path)) continue;
        current_file = directory::read_open(path);
        // Check all common in one file
        for (auto &value: common) {
            max_size = find_inside_file((*current_file), value.first);
            if (max_size > value.second.inputB.size()) {
                sub = value.first.substr(0, max_size);
                value.second.inputB = sub;
            }

            if (current_file->eof()) current_file->clear();
            current_file->seekg(0, ios::beg);
        }
        directory::read_close(current_file);
    }
}

void contig_diff::generate_output(const fs::path &directory_path, const map<std::string, contig_diff::Common> &common) {
    string outputPath = directory_path.string().append("/output.txt");
    ofstream *output = directory::write_open(outputPath, ios::trunc);
    (*output) << "Filename\tContig name\tContig value\tFind in B\tA - B\n";

    for (const auto &value: common) {
        (*output) << value.second.filename << "\t" << value.second.contigName << "\t" << value.first << "\t"
               << value.second.inputB << "\t" << value.first.substr(value.second.inputB.size()) << "\t\n";
    }

    directory::write_close(output);
}


