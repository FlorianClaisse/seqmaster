//
// Created by Florian Claisse on 02/05/2023.
//

#ifndef CONTIGDIFF_FASTADECODER_H
#define CONTIGDIFF_FASTADECODER_H

#include <filesystem>
#include <tuple>
#include <map>
#include <vector>
#include <functional>

namespace fasta {
    /** Transforme un fichier fasta vers un nouveau fichier en format fastaline. */
    int to_fasta_line(const std::filesystem::path &filePath);

    std::map<std::string, std::string> decode_fastaline(const std::filesystem::path &filePath);

    /** Permet de savoir si un fichier est de type fasta. */
    bool is_fasta_file(const std::filesystem::path &filePath);
    /** Permet de savoir si un fichier est de type result. */
    bool is_result_file(const std::filesystem::path &filePath);
    /** Permet de savoir si un fichier est de type fastaline. */
    bool is_fastaline_file(const std::filesystem::path &filePath);

    /** Dans un fichier de type fastaline permet de dire si un contig est présent. */
    bool find_contig(const std::filesystem::path &filePath, const std::string &contig);
    /** Dans un fichier de type fastaline permet de dire si tous les contigs sont présent ou non. */
    void find_contig(const std::filesystem::path &file_path, const std::map<std::string, std::string> &contigs, bool nucleic, std::function<void(const std::string&, const std::string&, const std::string&)> func);
    /** Dans un fichier de type fastaline permet de dire si tous les sont présent ou non avec un certains pourcentage d'erreur. */
    void find_contigs(const std::filesystem::path &file_path, const std::map<std::string, std::string> &contigs, int maxErrorPercentage, bool nucleic, std::function<void(const std::string&, const std::string&, const std::string&, double)> func);
}

#endif //CONTIGDIFF_FASTADECODER_H
