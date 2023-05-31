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

    int count_fasta_file(const std::filesystem::path &directoryPath);
    std::pair<std::filesystem::path, std::filesystem::path> two_first_fasta(const std::filesystem::path &directoryPath);

    int to_fastaline(const std::filesystem::path &directoryPath);
}

#endif //CONTIGDIFF_DIRECTORY_H
