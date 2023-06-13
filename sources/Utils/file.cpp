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

bool file::is_fastaline(const std::filesystem::path &path) {

    return have_extension(path, ".fastaline");
}

fs::path file::to_fastaline(const std::filesystem::path &filePath) {
    ifstream inputFile;
    inputFile.open(filePath);

    fs::path newPath{fs::path(filePath).replace_extension("fastaline")};
    ofstream outputFile(newPath, ios::trunc);

    string lineRead;
    bool first(true);
    while(getline(inputFile, lineRead)) {
        if (lineRead.empty()) continue;
        lineRead.erase(remove_if(lineRead.begin(), lineRead.end(), [](char c) { return c == '\n' || c == '\r'; }), lineRead.end());
        if (lineRead.at(0) == '>') {
            if (!first) outputFile << endl;
            outputFile << lineRead << endl;
            if (first) first = false;
        } else {
            transform(lineRead.begin(), lineRead.end(), lineRead.begin(), ::toupper);
            outputFile << lineRead;
        }
    }

    inputFile.close();
    outputFile.close();

    return newPath;
}