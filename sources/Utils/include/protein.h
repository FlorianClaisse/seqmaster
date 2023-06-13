//
// Created by Florian Claisse on 22/05/2023.
//

#ifndef CONTIG_PROTEIN_HPP
#define CONTIG_PROTEIN_HPP

#include <map>
#include <set>
#include <vector>
#include <string>

namespace protein {
    extern const std::map<char, std::vector<std::string>> amino_acide_codon_ref;
    extern const std::map<std::string, char> codon_ref_amino_acide;
    extern const std::set<std::string> codon_list;
    double rscu(const std::string &codon, const std::map<std::string, int> &total);
}

#endif //CONTIG_PROTEIN_HPP
