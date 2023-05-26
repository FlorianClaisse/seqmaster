//
// Created by Florian Claisse on 26/05/2023.
//

#ifndef CONTIG_PROTEIN_H
#define CONTIG_PROTEIN_H

#include <filesystem>

namespace find_all::protein {
    void search(const std::filesystem::path &inputA, const std::filesystem::path &inputB, const std::filesystem::path &output, int accept, int threads);
}

#endif //CONTIG_PROTEIN_H
