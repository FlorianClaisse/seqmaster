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
        {'G', {"GGT", "GGC", "GGA"}}
};