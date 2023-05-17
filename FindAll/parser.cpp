//
// Created by Florian Claisse on 15/05/2023.
//

#include <iostream>
#include <charconv>

#include "include/parser.h"
#include "include/find_all.h"

#include "../Utils/include/directory.h"

#define INPUTA "--inputA"
#define INPUTB "--inputB"
#define OUTPUT "--output"
#define TYPE "--type"
#define ACCEPT "--accept"
#define THREADS "--threads"

#define NUCL "nucl"
#define PROT "prot"

using namespace std;


namespace find_all {
    int usage();
    int check_options(const find_all::option &options);
}

// --inputA <path> --inputB <path> --type <nucl/prot> --output <path> [--accept <percentage>] [--threads <number>]
int find_all::parse(const std::vector<std::string_view> &argv) {
    if (argv.size() >= 8 && (argv.size() % 2) == 0) {
        if (argv[0] == INPUTA && argv[2] == INPUTB && argv[4] == TYPE && argv[6] == OUTPUT) {
            bool accept, thread;
            int acceptValue(100), threadValue(4);
            if (argv.size() >= 10 && argv.size() <= 12) {
                if (argv[8] == ACCEPT) {
                    accept = true;
                    auto result = from_chars(argv[9].data(), argv[9].data() + argv[9].size(), acceptValue);
                    if (result.ec == errc::invalid_argument) return usage();
                } else if (argv[8] == THREADS) {
                    thread = true;
                    auto result = from_chars(argv[9].data(), argv[9].data() + argv[9].size(), threadValue);
                    if (result.ec == errc::invalid_argument) return usage();
                } else return usage();
            }
            if (argv.size() == 12) {
                if (argv[10] == ACCEPT && !accept) {
                    auto result = from_chars(argv[11].data(), argv[11].data() + argv[11].size(), acceptValue);
                    if (result.ec == errc::invalid_argument) return usage();
                } else if (argv[10] == THREADS && !thread) {
                    auto result = from_chars(argv[11].data(), argv[11].data() + argv[11].size(), threadValue);
                    if (result.ec == errc::invalid_argument) return usage();
                } else return usage();
            }
            if (argv.size() > 12) return usage();

            if (string(argv[5]) != NUCL && string(argv[5]) != PROT) {
                cout << "Le type de fichier n'est pas valide.\n";
                return usage();
            }

            find_all::option options = {string(argv[1]), string(argv[3]), string(argv[7]), acceptValue, string(argv[5]) == NUCL, threadValue};
            if (check_options(options) == EXIT_FAILURE) return EXIT_FAILURE;

            return find_all::main(options);
        } else return usage();
    } else return usage();
}

int find_all::usage() {
    cout << "Find All" << endl
         << "Usage :" << endl
         << "./Contig " << "--findall" << " " << INPUTA << " <path> " << INPUTB << " <path> " << TYPE << " nucl/prot " << OUTPUT << " <path> " << "[" << ACCEPT << " <percentage>" << "]" << " [" << THREADS << " <value>]\n\n"
         << "\t" << INPUTA << "\tChemin vers le fichiers qui contient les contigs à trouver.\n\n"
         << "\t" << INPUTB << "\tChemin vers le dossier qui contient les fichiers ou il faut trouver les contigs.\n\n"
         << "\t" << TYPE << "\t\tLe type de fichier (nucl/prot).\n\n"
         << "\t" << OUTPUT << "\tChemin vers le dossier qui va contenir le(s) fichier(s) de sortie.\n\n"
         << "\t" << ACCEPT << "\tPermet de spécifier le pourcentage minimum pour accepter un contig comme reconnu.\n\n"
         << "\t" << THREADS << "\tDéfinit le nombre de threads que le programme peut utiliser\n\t\t\t(Default = 4)\n";
    return EXIT_SUCCESS;
}

int find_all::check_options(const find_all::option &options) {
    if (!directory::is_fasta_file(options.inputA)) {
        cout << "Path : " << options.inputA << " n'est pas un fichier fasta.\n";
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
    return EXIT_SUCCESS;
}