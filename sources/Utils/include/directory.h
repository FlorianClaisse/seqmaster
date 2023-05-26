//
// Created by Florian Claisse on 02/05/2023.
//

#ifndef CONTIGDIFF_DIRECTORY_H
#define CONTIGDIFF_DIRECTORY_H

#include <string>
#include <fstream>
#include <filesystem>
#include <ios>
#include <vector>

namespace directory {
    int create_directories(const std::filesystem::path &path);

    std::ifstream* read_open(const std::filesystem::path &filePath);
    std::ofstream* write_open(const std::filesystem::path &filePath, std::ios_base::openmode mode);

    void read_close(std::ifstream *file);
    void write_close(std::ofstream *file);

    bool have_extension(const std::filesystem::path &filePath, const std::vector<std::string> &extensions);

    bool is_fasta_file(const std::filesystem::path &filePath);
    bool is_result_file(const std::filesystem::path &filePath);
    bool is_fastaline_file(const std::filesystem::path &filePath);

    int count_fasta_file(const std::filesystem::path &directoryPath);
    std::pair<std::filesystem::path, std::filesystem::path> two_first_fasta(const std::filesystem::path &directoryPath);
}

#endif //CONTIGDIFF_DIRECTORY_H
