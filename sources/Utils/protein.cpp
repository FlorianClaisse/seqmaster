//
// Created by Florian Claisse on 22/05/2023.
//

#include "include/protein.h"

using namespace std;

const unordered_map<char, vector<string>> protein::amino_acide_codon_ref = {
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
        {'R', {"CGT, CGC", "CGA", "CGG", "AGA", "AGG"}},
        {'I', {"ATT", "ATC", "ATA"}},
        {'M', {"ATG"}},
        {'T', {"ACT", "ACC", "ACA", "ACG"}},
        {'N', {"AAT", "AAC"}},
        {'K', {"AAA", "AAG"}},
        {'V', {"GTT", "GTC", "GTA", "GTG"}},
        {'A', {"GCT", "GCC", "GCA", "GCG"}},
        {'D', {"GAT", "GAC"}},
        {'E', {"GAA", "GAG"}},
        {'G', {"GGT", "GGC", "GGA", "GGG"}},
};

const unordered_map<string, char> protein::codon_ref_amino_acide = {
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
        /*{'T', {"ACT", "ACC", "ACA", "ACG"}},
        {'N', {"AAT", "AAC"}},
        {'K', {"AAA", "AAG"}},
        {'V', {"GTT", "GTC", "GTA", "GTG"}},
        {'A', {"GCT", "GCC", "GCA", "GCG"}},
        {'D', {"GAT", "GAC"}},
        {'E', {"GAA", "GAG"}},
        {'G', {"GGT", "GGC", "GGA", "GGG"}},*/
};