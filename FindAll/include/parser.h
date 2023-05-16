//
// Created by Florian Claisse on 15/05/2023.
//

#ifndef CONTIG_FIND_ALL_PARSER_H
#define CONTIG_FIND_ALL_PARSER_H

#include <vector>
#include <filesystem>
#include <string_view>

namespace find_all {
    int parse(const std::vector<std::string_view> &argv);

    typedef struct {
        std::filesystem::path inputA;
        std::filesystem::path inputB;
        std::filesystem::path output;
        int accept;
        bool nucl; // nucl / prot
        int threads;
    } option;
}

#endif //CONTIG_FIND_ALL_PARSER_H
