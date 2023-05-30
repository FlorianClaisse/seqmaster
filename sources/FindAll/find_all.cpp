//
// Created by Florian Claisse on 02/05/2023.
//

#include <iostream>

#include "include/find_all.h"
#include "include/protein.h"
#include "include/nucleic.h"

#include "../Utils/include/termcolor.hpp"
#include "../Utils/include/file.h"
#include "../Utils/include/directory.h"

#define NUCL "nucl"
#define PROT "prot"

using namespace std;
namespace fs = std::filesystem;

int error_message(const string &message) {
    cout << termcolor::red << termcolor::bold
         << message
         << termcolor::reset;
    return -1;
}

int check_options(const fs::path &inputA, const fs::path &inputB, const fs::path &output, const string &type, int accept, int threads) {
    if (threads <= 0) return error_message("You want to run the program on a number of threads less than or equal to 0, go try again ;).\n");
    if (type != PROT && type != NUCL) return error_message("For the type, the value must be \"prit\" or \"nucl\" did you read the docs ?.\n");
    if (accept < 0 || accept > 100) return error_message("You just gave an acceptance percentage lower than 0 or higher than 100. That's not very clever.\n");
    if (!file::is_fasta(inputA)) return error_message("You need to give fasta file in inputA. Did you really read the doc.\n");
    if (!is_directory(inputB)) return error_message("Input B is not a folder or does not exist. Try again the next one is the right one.\n");
    if(!directory::create_directories(output)) return error_message("Unable to find/create the output folder are you sure you have given a valid path.\n");

    return 0;
}

int find_all::main(const fs::path &inputA, const fs::path &inputB, const fs::path &output, const string &type, int accept, int threads) {

    if (check_options(inputA, inputB, output, type, accept, threads) != 0) return -1;

    if (type == NUCL)
        nucleic::search(inputA, inputB, output, accept, threads);
    else if (type == PROT)
        protein::search(inputA, inputB, output, accept, threads);
    else
        return error_message("You just gave an invalid type. That's not very clever.\n");

    return 0;
}