//
// Created by Florian Claisse on 02/05/2023.
//

#include <iostream>
#include <algorithm>

#include "include/file.h"
#include "include/directory.h"

using namespace std;
namespace fs = std::filesystem;

int directory::create_directories(const std::filesystem::path &path) {
    error_code err;
    if (!fs::create_directories(path, err)) {
        if (exists(path)) return true;
        return false;
    }
    return true;
}

int directory::count_file(const std::filesystem::path &dirPath) {
    int cpt{0};
    for (const auto &v: fs::directory_iterator{dirPath}) {
        if (fs::is_regular_file(v))
            cpt++;
    }

    return cpt;
}

int directory::count_fasta_file(const std::filesystem::path &directoryPath) {
    int cpt(0);
    for(const auto &path: fs::directory_iterator(directoryPath)) {
        if (file::is_fasta(path)) cpt++;
    }
    return cpt;
}

std::pair<fs::path, fs::path> directory::two_first_fasta(const std::filesystem::path &directoryPath) {
    bool find_one(false);
    fs::path first;
    for (const auto &path: fs::directory_iterator(directoryPath)) {
        if (file::is_fasta(path)) {
            if (!find_one) {
                find_one = true;
                first = path;
            } else return {first, path};
        }
    }

    return {"", ""}; // Ne doit jamais arriver.
}

void directory::to_fastaline(const std::filesystem::path &directoryPath) {
    for (const auto &currentFile : fs::directory_iterator(directoryPath)) {
        if (!file::is_fasta(currentFile)) continue;
        file::to_fastaline(currentFile);
    }
}
