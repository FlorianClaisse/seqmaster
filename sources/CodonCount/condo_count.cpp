//
// Created by Florian Claisse on 04/05/2023.
//

#include <map>
#include <iostream>
#include <fstream>

#include "include/condo_count.h"

#include "../Utils/include/directory.h"

using namespace std;
namespace fs = std::filesystem;

namespace codon_count {
    void init_codon(map<string, pair<int, char>> &codon);
    void reset_codon(map<string, pair<int, char>> &codon);
}

int codon_count::main(const std::filesystem::path &input, const std::filesystem::path &output) {

    map<string, pair<int, char>> codon, total_codon;
    init_codon(codon);
    init_codon(total_codon);

    directory::to_fastaline(input);

    for (const auto &currentFile: fs::directory_iterator(input)) {
        if (!file::is_fastaline(currentFile)) continue;

        string fileName = currentFile.path().stem();
        ifstream *inputFile = file::read_open(currentFile.path());

        ofstream *outputFile = file::write_open(output.string().append("/" + fileName + ".txt"), ios::trunc);
        (*outputFile) << "Contig Name\tCodon\tAmino acids\tNumber\tPercentage\n";

        string lineRead, prot_name;
        unsigned long total(0), all_total(0);
        while (getline((*inputFile), lineRead)) {
            if (lineRead.at(0) == '>') prot_name = lineRead;
            else {
                for (unsigned long i = 0; i < lineRead.size(); i += 3) {
                    string sub = lineRead.substr(i, 3);
                    if (codon.find(sub) != codon.end()) {
                        codon[sub].first++;
                        total_codon[sub].first++;
                        total++;
                        all_total++;
                    } else {
                        cout << "Contig : " << prot_name << ", codon : " << sub << " unknown\n";
                    }
                }

                for (const auto &key: codon) {
                    if (key.second.first != 0) {
                        (*outputFile) << prot_name << "\t" << key.first << "\t" << key.second.second << "\t" << key.second.first << "\t" << (((double) key.second.first / (double) total) * 100) << "%\n";
                    }
                }
                total = 0;
                reset_codon(codon);
            }
        }

        file::write_close(outputFile);
        file::read_close(inputFile);

        ofstream *total_output = file::write_open(output.string().append("/" + fileName + "-total.txt"), ios::trunc);
        (*total_output) << "Codon\tAmino acids\tNumber\tPercentage\n";

        for (const auto &key : total_codon) {
            if (key.second.first == 0) continue;
            (*total_output) << key.first << "\t" << key.second.second << "\t" << key.second.first << "\t" << (((double)key.second.first / (double) all_total) * 100) << "%\n";
        }

        file::write_close(total_output);
        reset_codon(total_codon);
    }

    return EXIT_SUCCESS;
}

void codon_count::init_codon(map<string, pair<int, char>> &codon) {
    codon["TTT"] = {0, 'F'}; codon["TCT"] = {0, 'S'}; codon["TAT"] = {0, 'Y'}; codon["TGT"] = {0, 'C'};
    codon["TTC"] = {0, 'F'}; codon["TCC"] = {0, 'S'}; codon["TAC"] = {0, 'Y'}; codon["TGC"] = {0, 'C'};
    codon["TTA"] = {0, 'L'}; codon["TCA"] = {0, 'S'}; codon["TAA"] = {0, '*'}; codon["TGA"] = {0, '*'};
    codon["TTG"] = {0, 'L'}; codon["TCG"] = {0, 'S'}; codon["TAG"] = {0, '*'}; codon["TGG"] = {0, 'W'};

    codon["CTT"] = {0, 'L'}; codon["CCT"] = {0, 'P'}; codon["CAT"] = {0, 'H'}; codon["CGT"] = {0, 'R'};
    codon["CTC"] = {0, 'L'}; codon["CCC"] = {0, 'P'}; codon["CAC"] = {0, 'H'}; codon["CGC"] = {0, 'R'};
    codon["CTA"] = {0, 'L'}; codon["CCA"] = {0, 'P'}; codon["CAA"] = {0, 'Q'}; codon["CGA"] = {0, 'R'};
    codon["CTG"] = {0, 'L'}; codon["CCG"] = {0, 'P'}; codon["CAG"] = {0, 'Q'}; codon["CGG"] = {0, 'R'};

    codon["ATT"] = {0, 'I'}; codon["ACT"] = {0, 'T'}; codon["AAT"] = {0, 'N'}; codon["AGT"] = {0, 'S'};
    codon["ATC"] = {0, 'I'}; codon["ACC"] = {0, 'T'}; codon["AAC"] = {0, 'N'}; codon["AGC"] = {0, 'S'};
    codon["ATA"] = {0, 'I'}; codon["ACA"] = {0, 'T'}; codon["AAA"] = {0, 'K'}; codon["AGA"] = {0, 'R'};
    codon["ATG"] = {0, 'M'}; codon["ACG"] = {0, 'T'}; codon["AAG"] = {0, 'K'}; codon["AGG"] = {0, 'R'};

    codon["GTT"] = {0, 'V'}; codon["GCT"] = {0, 'A'}; codon["GAT"] = {0, 'D'}; codon["GGT"] = {0, 'G'};
    codon["GTC"] = {0, 'V'}; codon["GCC"] = {0, 'A'}; codon["GAC"] = {0, 'D'}; codon["GGC"] = {0, 'G'};
    codon["GTA"] = {0, 'V'}; codon["GCA"] = {0, 'A'}; codon["GAA"] = {0, 'E'}; codon["GGA"] = {0, 'G'};
    codon["GTG"] = {0, 'V'}; codon["GCG"] = {0, 'A'}; codon["GAG"] = {0, 'E'}; codon["GGG"] = {0, 'G'};
}

void codon_count::reset_codon(map<std::string, pair<int, char>> &codon) {
    for (auto &value: codon) {
        value.second.first = 0;
    }
}