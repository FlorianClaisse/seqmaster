//
// Created by Florian Claisse on 22/05/2023.
//

#ifndef CONTIG_PROTEIN_H
#define CONTIG_PROTEIN_H

#include <unordered_map>
#include <vector>
#include <string>

namespace protein {
    extern const std::unordered_map<char, std::vector<std::string>> amino_acide_codon_ref;
    extern const std::vector<std::string> codon_list;
}

#endif //CONTIG_PROTEIN_H
