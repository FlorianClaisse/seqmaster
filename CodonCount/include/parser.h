//
// Created by Florian Claisse on 15/05/2023.
//

#ifndef CONTIG_CODON_COUNT_PARSER_H
#define CONTIG_CODON_COUNT_PARSER_H

#include <filesystem>
#include <vector>
#include <string_view>

namespace codon_count {

    int parse(const std::vector<std::string_view> &argv);
    
    typedef struct {
        std::filesystem::path input;
        std::filesystem::path output;
    } option;
    
}

#endif //CONTIG_CODON_COUNT_PARSER_H