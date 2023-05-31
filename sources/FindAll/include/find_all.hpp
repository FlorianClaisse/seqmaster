//
// Created by Florian Claisse on 02/05/2023.
//

#ifndef CONTIG_FINDALL_FIND_ALL_HPP
#define CONTIG_FINDALL_FIND_ALL_HPP

#include <filesystem>

namespace find_all {
    int main(const std::filesystem::path &inputA, const std::filesystem::path &inputB, const std::filesystem::path &output,
             const std::string &type, int accept, int threads);

    struct param {
        std::filesystem::path inputA;
        std::filesystem::path inputB;
        std::filesystem::path output;
        bool nucl;
        int error_rate;
        int threads;
    };
}

#endif //CONTIG_FINDALL_FIND_ALL_HPP
