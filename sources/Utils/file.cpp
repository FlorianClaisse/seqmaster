//
// Created by Florian Claisse on 25/05/2023.
//

#include <seqan3/io/sequence_file/all.hpp>

#include "include/file.h"

using namespace std;
namespace fs = std::filesystem;

#define assert_message(e, m)                                    \
    do                                                          \
    {                                                           \
        if (!(e))                                               \
        {                                                       \
            std::cout << "Erreur : " << (m) << std::endl;       \
            std::abort();                                       \
        }                                                       \
    } while(0)

ifstream* file::read_open(const std::filesystem::path &filePath) {
    auto *file = new ifstream(filePath);
    assert_message(file->is_open(), "Impossible d'ouvrir le fichier");
    return file;
}

ofstream* file::write_open(const std::filesystem::path &filePath, ios_base::openmode mode) {
    auto *file = new ofstream(filePath, mode);
    assert_message(file->is_open(), "Impossible d'ouvrir le fichier");
    return file;
}

void file::read_close(ifstream *file) {
    file->close();
    delete file;
}
void file::write_close(ofstream *file) {
    file->close();
    delete file;
}

bool file::have_extension(const std::filesystem::path &path, const std::string &ext) {
    return is_regular_file(path) && path.extension() == ext;
}

bool file::have_extension(const std::filesystem::path &path, const std::vector<std::string> &exts) {
    return fs::is_regular_file(path) && find(exts.begin(), exts.end(), path.extension().string().erase(0, 1)) != exts.end();
}

bool file::is_fasta(const std::filesystem::path &path) {
    return have_extension(path, seqan3::format_fasta::file_extensions);
}