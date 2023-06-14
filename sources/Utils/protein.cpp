//
// Created by Florian Claisse on 22/05/2023.
//

#include "include/protein.h"

using namespace std;

double protein::rscu(const std::string &codon, const std::map<std::string, int> &total) {
    vector<string> codons{amino_acide_codon_ref.find(codon_ref_amino_acide.find(codon)->second)->second};
    size_t n_i = codons.size();
    int x_i_j = total.find(codon)->second;

    size_t frac = n_i * x_i_j;

    int sum{0};
    for (int i = 0; i < n_i; i++)
        sum += total.find(codons[i])->second;

    return (double) frac / (double)sum;
}

const map<char, vector<string>> protein::amino_acide_codon_ref = {
        {'F', {"TTT", "TTC"}},
        {'L', {"TTA", "TTG", "CTT", "CTC", "CTA", "CTG"}},
        {'S', {"TCT", "TCC", "TCA", "TCG", "AGT", "AGC"}},
        {'Y', {"TAT", "TAC"}},
        {'*', {"TAA", "TAG", "TGA"}},
        {'C', {"TGT", "TGC"}},
        {'W', {"TGG"}},
        {'P', {"CCT", "CCC", "CCA", "CCG"}},
        {'H', {"CAT", "CAC"}},
        {'Q', {"CAA", "CAG"}},
        {'R', {"CGT", "CGC", "CGA", "CGG", "AGA", "AGG"}},
        {'I', {"ATT", "ATC", "ATA"}},
        {'M', {"ATG"}},
        {'T', {"ACT", "ACC", "ACA", "ACG"}},
        {'N', {"AAT", "AAC"}},
        {'K', {"AAA", "AAG"}},
        {'V', {"GTT", "GTC", "GTA", "GTG"}},
        {'A', {"GCT", "GCC", "GCA", "GCG"}},
        {'D', {"GAT", "GAC"}},
        {'E', {"GAA", "GAG"}},
        {'G', {"GGT", "GGC", "GGA", "GGG"}}
};

const map<string, char> protein::codon_ref_amino_acide = {
        {"TTT", 'F'}, {"TTC", 'F'},
        {"TTA", 'L'}, {"TTG", 'L'}, {"CTT", 'L'}, {"CTC", 'L'}, {"CTA", 'L'}, {"CTG", 'L'},
        {"TCT", 'S'}, {"TCC", 'S'}, {"TCA", 'S'}, {"TCG", 'S'}, {"AGT", 'S'}, {"AGC", 'S'},
        {"TAT", 'Y'}, {"TAC", 'Y'},
        {"TAA", '*'}, {"TAG", '*'}, {"TGA", '*'},
        {"TGT", 'C'}, {"TGC", 'C'},
        {"TGG", 'W'},
        {"CCT", 'P'}, {"CCC", 'P'}, {"CCA", 'P'}, {"CCG", 'P'},
        {"CAT", 'H'}, {"CAC", 'H'},
        {"CAA", 'Q'}, {"CAG", 'Q'},
        {"CGT", 'R'}, {"CGC", 'R'}, {"CGA", 'R'}, {"CGG", 'R'}, {"AGA", 'R'}, {"AGG", 'R'},
        {"ATT", 'I'}, {"ATC", 'I'}, {"ATA", 'I'},
        {"ATG", 'W'},
        {"ACT", 'T'}, {"ACC", 'T'}, {"ACA", 'T'}, {"ACG", 'T'},
        {"AAT", 'N'}, {"AAC", 'N'},
        {"AAA", 'K'}, {"AAG", 'K'},
        {"GTT", 'V'}, {"GTC", 'V'}, {"GTA", 'V'}, {"GTG", 'V'},
        {"GCT", 'A'}, {"GCC", 'A'}, {"GCA", 'A'}, {"GCG", 'A'},
        {"GAT", 'D'}, {"GAC", 'D'},
        {"GAA", 'E'}, {"GAG", 'E'},
        {"GGT", 'G'}, {"GGC", 'G'}, {"GGA", 'G'}, {"GGG", 'G'}
};

const std::set<std::string> protein::codon_list = {
        "TTT", "TTC",
        "TTA", "TTG", "CTT", "CTC", "CTA", "CTG",
        "TCT", "TCC", "TCA", "TCG", "AGT", "AGC",
        "TAT", "TAC",
        "TAA", "TAG", "TGA",
        "TGT", "TGC",
        "TGG",
        "CCT", "CCC", "CCA", "CCG",
        "CAT", "CAC",
        "CAA", "CAG",
        "CGT", "CGC", "CGA", "CGG", "AGA", "AGG",
        "ATT", "ATC", "ATA",
        "ATG",
        "ACT", "ACC", "ACA", "ACG",
        "AAT", "AAC",
        "AAA", "AAG",
        "GTT", "GTC", "GTA", "GTG",
        "GCT", "GCC", "GCA", "GCG",
        "GAT", "GAC",
        "GAA", "GAG",
        "GGT", "GGC", "GGA", "GGG"
};