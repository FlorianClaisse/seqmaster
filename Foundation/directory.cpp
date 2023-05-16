//
// Created by Florian Claisse on 02/05/2023.
//

#include <iostream>

#include "include/directory.h"

#define assert_message(e, m)                                    \
    do                                                          \
    {                                                           \
        if (!(e))                                               \
        {                                                       \
            std::cout << "Erreur : " << (m) << std::endl;       \
            std::abort();                                       \
        }                                                       \
    } while(0)

using namespace std;
namespace fs = std::filesystem;

std::ifstream* directory::read_open(const std::filesystem::path &filePath) {
    auto *file = new ifstream(filePath);
    assert_message(file->is_open(), "Erreur : Impossible d'ouvrir le fichier");
    return file;
}

std::ofstream* directory::write_open(const std::filesystem::path &filePath, ios_base::openmode mode) {
    auto *file = new ofstream(filePath, mode);
    assert_message(file->is_open(), "Erreur : Impossible d'ouvrir le fichier");
    return file;
}

void directory::read_close(std::ifstream *file) {
    file->close();
    delete file;
}
void directory::write_close(std::ofstream *file) {
    file->close();
    delete file;
}

bool directory::have_extension(const fs::path &filePath, const vector<string> &extensions) {
    return fs::is_regular_file(filePath) && find(extensions.begin(), extensions.end(), filePath.extension()) != extensions.end();
}

bool directory::is_fasta_file(const fs::path &filePath) {
    vector<string> extensions = {".fasta", ".fna", ".faa", ".ffn", ".fa", ".fas"};
    return have_extension(filePath, extensions);
}

bool directory::is_result_file(const std::filesystem::path &filePath) {
    return fs::is_regular_file(filePath) && (filePath.string().find("-result.fasta") != string::npos);
}

bool directory::is_fastaline_file(const fs::path &filePath) {
    return fs::is_regular_file(filePath) && filePath.extension() == ".fastaline";
}

int directory::count_fasta_file(const std::filesystem::path &directoryPath) {
    int cpt(0);
    for(const auto &path: fs::directory_iterator(directoryPath)) {
        if (is_fastaline_file(path)) cpt++;
    }
    return cpt;
}

std::pair<fs::path, fs::path> directory::two_first_fasta(const std::filesystem::path &directoryPath) {
    bool find_one(false);
    fs::path first;
    for (const auto &path: fs::directory_iterator(directoryPath)) {
        if (is_fastaline_file(path)) {
            if (!find_one) {
                find_one = true;
                first = path;
            } else return {first, path};
        }
    }

    return {"", ""};
}
