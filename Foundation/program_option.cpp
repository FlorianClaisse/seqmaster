//
// Created by Florian Claisse on 02/05/2023.
//

#include <filesystem>
#include <charconv>

#include "include/program_option.h"
#include "../FindAll/find_all.h"
#include "../CodonCount/condo_count.h"
#include "../ContigDiff/contig_diff.h"

using namespace std;
namespace fs = std::filesystem;

int pathDosentExists(const string_view &path) {
    cout << "Path : " << path << " n'existe pas ou n'est pas accessible.\n" ;
    return EXIT_FAILURE;
}

// ./Contig [program_name] ...
int program_option::parse(int argc, char **argv) {
    vector<string_view> args(argv, argv + argc);
    if (args.size() < 2) {
        cout << "La ligne de commande doit au moins contenir le nom du program à éxécuter." << endl;
        return usage();
    }

    vector<string_view> sub_args(argv + 2, argv + argc);
    if (args[1] == FINDALL) {
        return parse_find_all(sub_args);
    } else if (args[1] == CODONCOUNT) {
        return parse_codon_count(sub_args);
    } else if (args[1] == CONTIGDIFF) {
        return parse_contig_diff(sub_args);
    }

    return usage();
}

int program_option::usage() {
    cout << "Usage :" << endl
    << "./Contig [program_name] ..." << endl
    << "\t" << FINDALL << "\tProgramme qui permet à partir d'un fichier A (type nucl ou prot), de trouver si il sont présent dans tous les fichiers du dossier B.\n"
    << "\t" << CODONCOUNT << "\tProgramme qui permet à partir d'un fichier d'entrée A de compter le nombre de chaque codon pour chaque contig.\n";

    return EXIT_SUCCESS;
}

// --inputA <path> --inputB <path> --type <nucl/prot> --output <path> [--accept <percentage>] [--threads <value>]
int program_option::parse_find_all(const vector<string_view> &argv) {
    if (argv.size() < 8 || (argv.size() % 2) != 0) return find_all_usage();
    if (argv[0] != INPUTA || argv[2] != INPUTB || argv[4] != TYPE || argv[6] != OUTPUT) return find_all_usage();

    bool accept, thread;
    int acceptValue(100), threads(4);
    if (argv.size() >= 10 && argv.size() <= 12) {
        if (argv[8] == ACCEPT) {
            accept = true;
            auto result = from_chars(argv[9].data(), argv[9].data() + argv[9].size(), acceptValue);
            if (result.ec == errc::invalid_argument) return find_all_usage();
        } else if (argv[8] == THREADS) {
            thread = true;
            auto result = from_chars(argv[9].data(), argv[9].data() + argv[9].size(), threads);
            if (result.ec == errc::invalid_argument) return find_all_usage();
        } else return find_all_usage();
    }
    if (argv.size() == 12) {
        if (argv[10] == ACCEPT && !accept) {
            auto result = from_chars(argv[11].data(), argv[11].data() + argv[11].size(), acceptValue);
            if (result.ec == errc::invalid_argument) return find_all_usage();
        } else if (argv[10] == THREADS && !thread) {
            auto result = from_chars(argv[11].data(), argv[11].data() + argv[11].size(), threads);
            if (result.ec == errc::invalid_argument) return find_all_usage();
        } else return find_all_usage();
    }
    if (argv.size() > 12) return find_all_usage();

    if (!fs::exists(argv[1])) return pathDosentExists(argv[1]);

    if (!fs::exists(argv[3])) return pathDosentExists(argv[3]);

    if (string(argv[5]) != NUCLEIC && string(argv[5]) != PROTEIN) {
        cout << "Le type de fichier n'est pas valide." << endl;
        return EXIT_FAILURE;
    }

    if (!fs::exists(argv[7])) return pathDosentExists(argv[7]);

    FindAll options = {string(argv[1]), string(argv[3]), string(argv[7]), acceptValue, string(argv[5]) == NUCLEIC, threads};
    return find_all::start(options);
}


int program_option::find_all_usage() {
    cout << "Find All" << endl
    << "Usage :" << endl
    << "./Contig " << FINDALL << " " << INPUTA << " <path> " << INPUTB << " <path> " << TYPE << " nucl/prot " << OUTPUT << " <path> " << "[" << ACCEPT << " <percentage>" << "]" << " [" << THREADS << " <value>] \n\n"
    << "\t" << INPUTA << "\tChemin vers le fichiers qui contient les contigs à trouver.\n\n"
    << "\t" << INPUTB << "\tChemin vers le dossier qui contient les fichiers ou il faut trouver les contigs.\n\n"
    << "\t" << TYPE << "\t\tLe type de fichier (nucl/prot).\n\n"
    << "\t" << OUTPUT << "\tChemin vers le dossier qui va contenir le(s) fichier(s) de sortie.\n\n"
    << "\t" << ACCEPT << "\tPermet de spécifier le pourcentage minimum pour accepter un contig comme reconnu.\n\n"
    << "\t" << THREADS << "\tDéfinit le nombre de threads que le programme peut utiliser\n\t\t\t(Default = 4)\n";
    return EXIT_SUCCESS;
}

// --inputA <path> --output <path>
int program_option::parse_codon_count(const vector<string_view> &argv) {
    if (argv.size() != 4) return codon_count_usage();

    if(argv[0] != INPUTA || argv[2] != OUTPUT) return codon_count_usage();

    if (!fs::exists(argv[1])) return pathDosentExists(argv[1]);

    if (!fs::exists(argv[3])) return pathDosentExists(argv[3]);

    CodonCount options = {string(argv[1]), string(argv[3])};
    return codon_count::start(options);
}

int program_option::codon_count_usage() {
    cout << "Codon Count" << endl
    << "Usage :" << endl
    << "./Contig " << CODONCOUNT << " " << INPUTA << " <path> " << OUTPUT << " <path>" << endl
    << "\t" << INPUTA << "\tChemin vers le fichier qui contient les contigs ou il faut compter les codons.\n\n"
    << "\t" << OUTPUT << "\tChemin vers le dossier de sortie.\n";
    return EXIT_SUCCESS;
}

// --inputA <path> --inputB <path> --output <path> --accept <percentage> [--threads <number>]
int program_option::parse_contig_diff(const std::vector<std::string_view> &argv) {
    if (argv.size() < 8 || (argv.size() % 2) != 0) return contig_diff_usage();
    if (argv[0] != INPUTA || argv[2] != INPUTB || argv[4] != OUTPUT || argv[6] != ACCEPT) return contig_diff_usage();

    int threads(4);
    if (argv.size() == 10) {
        if (argv[8] == THREADS) {
            auto result = from_chars(argv[9].data(), argv[9].data() + argv[9].size(), threads);
            if (result.ec == errc::invalid_argument) return contig_diff_usage();
        } else return contig_diff_usage();
    } else if (argv.size() > 10) return contig_diff_usage();

    if (!fs::exists(argv[1])) return pathDosentExists(argv[1]);

    if (!fs::exists(argv[3])) return pathDosentExists(argv[3]);

    if (!fs::exists(argv[5])) return pathDosentExists(argv[5]);

    int accept(100);
    auto result = from_chars(argv[7].data(), argv[7].data() + argv[7].size(), accept);
    if (result.ec == errc::invalid_argument) return contig_diff_usage();

    program_option::ContigDiff options = {string(argv[1]), string(argv[3]), string(argv[5]), accept, threads};
    return contig_diff::start(options);
}

int program_option::contig_diff_usage() {
    cout << "Contig Diff" << endl
    << "Usage :" << endl
    << "./Contig " << CONTIGDIFF << " " << INPUTA << " <path> " << INPUTB << " <path> " << OUTPUT << " <path> " << ACCEPT << " <percentage> " << "[" << THREADS << " <number>]\n"
    << "\t" << INPUTA << "\tChemin vers le dossier A.\n\n"
    << "\t" << INPUTB << "\tChemin vers le dossier B.\n\n"
    << "\t" << OUTPUT << "\tChemin vers le dossier qui va contenir les fichiers de sortie.\n\n"
    << "\t" << ACCEPT << "\tLe pourcentage minimum de similitude accepté.\n\n"
    << "\t" << THREADS << "\tNombre de threads que le programme va utiliser (default=4)\n";
    return EXIT_SUCCESS;
}



