//
// Created by Florian Claisse on 02/05/2023.
//

#ifndef CONTIGDIFF_PROGRAM_OPTION_H
#define CONTIGDIFF_PROGRAM_OPTION_H

#include <iostream>
#include <vector>
#include <filesystem>

// Program name
#define FINDALL "--findAll"

// Commande option
#define INPUTA "--inputA"
#define INPUTB "--inputB"
#define OUTPUT "--output"
#define ACCEPT "--accept"

namespace program_option {
    /** Permet à partir d'une ligne de commande de savoir quel programme est demandé et de l'envoyé vers le bon parser. */
    int parse(int argc, char *argv[]);
    /** Parse la ligne de commande reconnu comme etant pour le programme find all. Une fois la commande persé correctement le programm est lancé. */
    int parse_find_all(const std::vector<std::string_view> &argv);

    /** Affiche les usage pour la ligne de commande des programmes. */
    int usage();
    /** Affiche les usage pour le programme find_all. */
    int find_all_usage();

    /** structure contenant les options necessaire pour le programme find all. */
    typedef struct {
        std::filesystem::path inputA;
        std::filesystem::path inputB;
        std::filesystem::path output;
        int accept;
    } FindAll;
}

#endif //CONTIGDIFF_PROGRAM_OPTION_H
