//
// Created by Florian Claisse on 16/05/2023.
//

#include <iostream>
#include <charconv>

#include "include/parser.h"
#include "include/contig_diff.h"

#define INPUTA "--inputA"
#define INPUTB "--inputB"
#define OUTPUT "--output"
#define ACCEPT "--accept"
#define THREADS "--threads"

using namespace std;

namespace contig_diff {
    int usage();
    int check_option(const contig_diff::option &options);
}

// --inputA <path> --inputB <path> --output <path> [--accept <percentage>] [--threads <value>]
int contig_diff::parse(std::vector<std::string_view> &argv) {
    if (argv.size() >= 6 && argv.size() <= 10 && (argv.size() % 2) == 0) {
        if (argv[0] == INPUTA && argv[2] == INPUTB && argv[4] == OUTPUT) {
            bool accept(false), thread(false);
            int acceptValue(100), threadValue(4);
            if (argv.size() >= 8) {
                if (argv[6] == THREADS) {
                    thread = true;
                    auto result = from_chars(argv[7].data(), argv[7].data() + argv[7].size(), threadValue);
                    if (result.ec == errc::invalid_argument) return usage();
                } else if (argv[6] == ACCEPT) {
                    accept = true;
                    auto result = from_chars(argv[7].data(), argv[7].data() + argv[7].size(), acceptValue);
                    if (result.ec == errc::invalid_argument) return usage();
                } else return usage();
            }
            if (argv.size() == 10) {
                if (argv[8] == THREADS && !thread) {
                    auto result = from_chars(argv[9].data(), argv[9].data() + argv[9].size(), threadValue);
                    if (result.ec == errc::invalid_argument) return usage();
                } else if (argv[8] == ACCEPT && !accept) {
                    auto result = from_chars(argv[9].data(), argv[9].data() + argv[9].size(), acceptValue);
                    if (result.ec == errc::invalid_argument) return usage();
                } else return usage();
            }

            contig_diff::option options = {string(argv[1]), string(argv[3]), string(argv[5]), acceptValue, acceptValue};
            if (check_option(options) == EXIT_FAILURE) return EXIT_FAILURE;

            //return contig_diff::main(options);
        } else return usage();
    } else return usage();
}

int contig_diff::check_option(const contig_diff::option &options) {
    if (!is_directory(options.inputA)) {
        cout << "Path : " << options.inputA << " n'est pas un dossier.\n";
        return EXIT_FAILURE;
    }
    if (!is_directory(options.inputB)) {
        cout << "Path : " << options.inputB << " n'est pas un dossier.\n";
        return EXIT_FAILURE;
    }
    if (!is_directory(options.output)) {
        cout << "Path : " << options.output << " n'est pas un dossier.\n";
        return EXIT_FAILURE;
    }

    if (options.accept < 0 || options.accept > 100) {
        cout << "Accept value must be between [0 - 100], current value is : " << options.accept << endl;
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}

int contig_diff::usage() {
    cout << "Contig Diff" << endl
         << "Usage :" << endl
         << "./Contig " << "--contigdiff" << " " << INPUTA << " <path> " << INPUTB << " <path> " << OUTPUT << " <path> " << "[" << THREADS << " <number>]\n"
         << "\t" << INPUTA << "\tChemin vers le dossier A.\n\n"
         << "\t" << INPUTB << "\tChemin vers le dossier B.\n\n"
         << "\t" << OUTPUT << "\tChemin vers le dossier qui va contenir les fichiers de sortie.\n\n"
         << "\t" << THREADS << "\tNombre de threads que le programme va utiliser (default=4)\n";

    return EXIT_SUCCESS;
}