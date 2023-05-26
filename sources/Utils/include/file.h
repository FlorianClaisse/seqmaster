//
// Created by Florian Claisse on 25/05/2023.
//

#ifndef CONTIG_FILE_H
#define CONTIG_FILE_H

#include <string>
#include <fstream>
#include <filesystem>
#include <vector>


namespace file {
    std::ofstream* write_open(const std::filesystem::path &path, std::ios_base::openmode mode);
    std::ifstream* read_open(const std::filesystem::path &path);

    void write_close(std::ofstream *path);
    void read_close(std::ifstream *path);

    bool have_extension(const std::filesystem::path &path, const std::string &ext);
    bool have_extension(const std::filesystem::path &path, const std::vector<std::string> &exts);
    bool is_fasta(const std::filesystem::path &path);
}

#endif //CONTIG_FILE_H
