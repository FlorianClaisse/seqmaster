//
// Created by Florian Claisse on 10/05/2023.
//

#ifndef CONTIG_CONTIG_DIFF_H
#define CONTIG_CONTIG_DIFF_H

#include <string>
#include <filesystem>
#include "parser.h"

namespace contig_diff {
    typedef struct {
        std::string filename;
        std::string contigName;
        std::string inputB;
        double percentage;
    } Common;

    int main(const std::filesystem::path &inputA, const std::filesystem::path &inputB, const std::filesystem::path &output,
             const std::string &type, int accept, int threads);
}

#endif //CONTIG_CONTIG_DIFF_H
