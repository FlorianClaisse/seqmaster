//
// Created by Florian Claisse on 16/05/2023.
//

#ifndef CONTIG_CONTIG_DIFF_PARSER_H
#define CONTIG_CONTIG_DIFF_PARSER_H

#include <filesystem>
#include <vector>
#include <string_view>

namespace contig_diff {

    int parse(std::vector<std::string_view> &argv);

    typedef struct {
        std::filesystem::path inputA;
        std::filesystem::path inputB;
        std::filesystem::path output;
        int accept;
        int threads;
    } option;
}

#endif //CONTIG_CONTIG_DIFF_PARSER_H
