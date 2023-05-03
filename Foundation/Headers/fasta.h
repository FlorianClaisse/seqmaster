//
// Created by Florian Claisse on 02/05/2023.
//

#ifndef CONTIGDIFF_FASTADECODER_H
#define CONTIGDIFF_FASTADECODER_H

#include <filesystem>
#include <tuple>
#include <map>

namespace fasta {
    /** Transforme un fichier fasta vers un nouveau fichier en format fastaline. */
    int to_fasta_line(const std::filesystem::path &filePath);

    /** Dans un fichier de type fastaline permet de dire si un contig est présent. */
    bool find_contig(const std::filesystem::path &filePath, const std::string &contig);
    /** Dans un fichier de type fastaline permet de dire si tous les contigs sont présent ou non. */
    std::map<std::string, bool>find_contigs(const std::filesystem::path &file_path,  const std::vector<std::string> &contigs_value);
}

#endif //CONTIGDIFF_FASTADECODER_H
