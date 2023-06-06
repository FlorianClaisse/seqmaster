//
// Created by Florian Claisse on 10/05/2023.
//

#ifndef CONTIG_CONTIG_DIFF_H
#define CONTIG_CONTIG_DIFF_H

#include <filesystem>

namespace contig_diff {

    int main(const std::filesystem::path &inputA, const std::filesystem::path &inputB, const std::filesystem::path &output,
             const std::string &type, int accept, int threads, unsigned long min_size);

    struct param {
        std::filesystem::path inputA;
        std::filesystem::path inputB;
        std::filesystem::path output;
        bool nucl;
        int error_rate;
        int threads;
        unsigned long min_size;
    };
}

#endif //CONTIG_CONTIG_DIFF_H
