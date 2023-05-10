//
// Created by Florian Claisse on 02/05/2023.
//

#ifndef CONTIGDIFF_PROGRAM_OPTION_H
#define CONTIGDIFF_PROGRAM_OPTION_H

#include <iostream>
#include <vector>
#include <filesystem>

// Program name
#define FINDALL "--findall"
#define CODONCOUNT "--codoncount"
#define CONTIGDIFF "--contigdiff"

// Commande option
#define THREADS "--threads"
#define INPUTA "--inputA"
#define INPUTB "--inputB"
#define TYPE "--type"
#define OUTPUT "--output"
#define ACCEPT "--accept"

#define PROTEIN "prot"
#define NUCLEIC "nucl"

namespace program_option {
    /** Permet à partir d'une ligne de commande de savoir quel programme est demandé et de l'envoyé vers le bon parser. */
    int parse(int argc, char *argv[]);
    /** Parse la ligne de commande reconnu comme etant pour le programme find all. Une fois la commande persé correctement le programm est lancé. */
    int parse_find_all(const std::vector<std::string_view> &argv);
    /** Parse le ligne de commande reconnu comme etant pour le programme codon count. Une fois la commande parsé correctement le program est lancé. */
    int parse_codon_count(const std::vector<std::string_view> &argv);
    /** Parse la ligne de commande reconnu comme etant cela pour le programme contig diff. Une fois la ligne de commande parsé, le programme est lancé. */
    int parse_contig_diff(const std::vector<std::string_view> &argv);

    /** Affiche les usage pour la ligne de commande des programmes. */
    int usage();
    /** Affiche les usage pour le programme find_all. */
    int find_all_usage();
    /** Affiche les usages pout le programme count_count */
    int codon_count_usage();
    /** Affiche les usages pour le programme contig diff. */
    int contig_diff_usage();

    /** structure contenant les options necessaire pour le programme find all. */
    typedef struct {
        std::filesystem::path inputA;
        std::filesystem::path inputB;
        std::filesystem::path output;
        int accept;
        bool nucl; /* nucl | prot */
        int nb_thread;
    } FindAll;

    typedef struct {
        std::filesystem::path inputA;
        std::filesystem::path output;
    } CodonCount;

    typedef struct {
        std::filesystem::path inputA;
        std::filesystem::path inputB;
        std::filesystem::path output;
        int accept;
        int threads;
    } ContigDiff;
}

#endif //CONTIGDIFF_PROGRAM_OPTION_H
