//
// Created by Florian Claisse on 26/05/2023.
//

#ifndef CONTIG_FINDALL_PROTEIN_HPP
#define CONTIG_FINDALL_PROTEIN_HPP

#include <filesystem>

namespace find_all::protein {
    void search(const std::filesystem::path &inputA, const std::filesystem::path &inputB, const std::filesystem::path &output, int accept, int threads);
}

#endif //CONTIG_FINDALL_PROTEIN_HPP
