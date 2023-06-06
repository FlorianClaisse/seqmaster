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

#include <seqan3/io/sequence_file/all.hpp>

#include "file.h"

namespace directory {
    int create_directories(const std::filesystem::path &path);

    int count_fasta_file(const std::filesystem::path &directoryPath);
    std::pair<std::filesystem::path, std::filesystem::path> two_first_fasta(const std::filesystem::path &directoryPath);

    int to_fastaline(const std::filesystem::path &directoryPath);

    template<typename traits_t, typename record_t>
    void decode_all_fasta(const std::filesystem::path &dir, std::vector<std::vector<record_t>> &all_records) {
        for (const auto &path: std::filesystem::directory_iterator(dir)) {
            if (!file::is_fasta(path)) continue;
            seqan3::sequence_file_input<traits_t> f_in{path};
            std::vector<record_t> records;
            std::ranges::copy(f_in, std::back_inserter(records));
            all_records.push_back(records);
        }
    }
}

#endif //CONTIGDIFF_DIRECTORY_H
