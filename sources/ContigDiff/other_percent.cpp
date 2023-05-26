//
// Created by Florian Claisse on 22/05/2023.
//

#include <iostream>
#include <filesystem>

#include "../Utils/include/directory.h"

#include "include/other_percent.h"

using namespace std;
namespace fs = std::filesystem;
/*
namespace contig_diff::other_percent {
    void search_with_2_file(const fs::path &first_path, const fs::path &second_path, map <string, Common> &common);
    unsigned long find_inside_file(ifstream &test_file, const string &word);
    void check_on_directory(const fs::path &directory_path, map<string, Common> &common);
    void verif_on_directory(const fs::path &directory_path, map<string, Common> &common);
}

map<string, contig_diff::Common> contig_diff::other_percent::main(const contig_diff::option &options) {
    cout << "Start finding all common inside A.\n";
    map<string, Common> common;

    pair<fs::path, fs::path> two_first = directory::two_first_fasta(options.inputA);
    search_with_2_file(two_first.first, two_first.second, common);

    for (const auto &path: fs::directory_iterator(options.inputA)) {
        if (path == two_first.first || !directory::is_fastaline_file(path)) continue;
        search_with_2_file(path, two_first.first, common);
    }

    return common;
}*/