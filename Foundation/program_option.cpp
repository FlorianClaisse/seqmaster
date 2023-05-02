//
// Created by Florian Claisse on 02/05/2023.
//

#include <filesystem>
#include <charconv>

#include "Headers/program_option.h"
#include "../FindAll/find_all.h"

// ./Contig [program_name] ...
int program_option::parse(int argc, char **argv) {
    std::vector<std::string_view> args(argv, argv + argc);
    if (args.size() < 2) {
        std::cout << "La ligne de commande doit au moins contenir le nom du program à éxécuter." << std::endl;
        return usage();
    }

    std::vector<std::string_view> sub_args(argv + 2, argv + argc);
    if (args[1] == FINDALL) {
        return parse_find_all(sub_args);
    }

    return usage();
}

int program_option::usage() {
    std::cout << "Usage :" << std::endl
    << "./Contig [program_name] ..." << std::endl
    << "\t" << FINDALL << "\tProgramme qui permet à partir d'un fichier A qui contient des contigs, de trouver si il sont présent dans tous les fichiers du dossier B."
    << std::endl;

    return EXIT_SUCCESS;
}

// --inputA <path> --inputB <path> [--output <path>] [--accept <percentage>]
int program_option::parse_find_all(const std::vector<std::string_view> &argv) {
    if (argv.size() < 4 || (argv.size() % 2) != 0) return find_all_usage();
    if (argv[0] != INPUTA || argv[2] != INPUTB) return find_all_usage();

    bool output(false), accept(false);
    std::string outputPath;
    int acceptValue(100);
    if (argv.size() >= 6) {
        if (argv[4] == OUTPUT) {
            output = true;
            outputPath = std::string(argv[5]);
        }
        else if (argv[4] == ACCEPT) {
            accept = true;
            auto result = std::from_chars(argv[5].data(), argv[5].data() + argv[5].size(), acceptValue);
            if (result.ec == std::errc::invalid_argument) return find_all_usage();
        }
        else return find_all_usage();
    }
    if (argv.size() == 8) {
        if (argv[6] == OUTPUT && !output) {
            outputPath = std::string(argv[5]);
        }
        else if (argv[6] == ACCEPT && !accept) {
            auto result = std::from_chars(argv[5].data(), argv[5].data() + argv[5].size(), acceptValue);
            if (result.ec == std::errc::invalid_argument) return find_all_usage();
        }
        else return find_all_usage();
    }
    if (argv.size() > 8) return find_all_usage();

    if (!std::filesystem::exists(argv[1])) { std::cout << "Le fichier d'entrée A n'existe pas ou n'est pas accessible." << std::endl; }
    if (!std::filesystem::exists(argv[3])) { std::cout << "Le dossier d'entrée B n'exsite pas ou n'est pas accessible." << std::endl; }

    FindAll options = {std::string(argv[1]), std::string(argv[3]), outputPath, acceptValue};
    return find_all::start(options);
}


int program_option::find_all_usage() {
    std::cout << "Find All" << std::endl
    << "Usage :" << std::endl
    << "\t" << INPUTA << "\tChemin vers le fichiers qui contient les contigs à trouver." << std::endl
    << "\t" << INPUTB << "\tChemin vers le dossier qui contient les fichiers ou il faut trouver les contigs." << std::endl
    << "\t" << OUTPUT << "\tChemin vers le dossier qui va contenir le/les fichier(s) de sortie." << std::endl
    << "\t" << ACCEPT << "\tPermet de spécifier le pourcentage minimum pour accepter un contig comme reconnu." << std::endl;
    return EXIT_SUCCESS;
}

