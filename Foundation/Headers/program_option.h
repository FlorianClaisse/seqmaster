//
// Created by Florian Claisse on 02/05/2023.
//

#ifndef CONTIGDIFF_PROGRAM_OPTION_H
#define CONTIGDIFF_PROGRAM_OPTION_H

#include <iostream>
#include <vector>

// Program name
#define FINDALL "--findAll"

// Commande option
#define INPUTA "--inputA"
#define INPUTB "--inputB"
#define OUTPUT "--output"
#define ACCEPT "--accept"

namespace program_option {
    int parse(int argc, char *argv[]);
    int parse_find_all(const std::vector<std::string_view> &argv);

    int usage();
    int find_all_usage();

    typedef struct {
        std::string inputA;
        std::string inputB;
        std::string output;
        int accept;
    } FindAll;
}

#endif //CONTIGDIFF_PROGRAM_OPTION_H
