//
// Created by Florian Claisse on 02/05/2023.
//

#include <filesystem>
#include <charconv>

#include "include/program_option.h"
#include "../FindAll/include/find_all.h"
#include "../CodonCount/condo_count.h"
#include "../ContigDiff/contig_diff.h"

using namespace std;
namespace fs = std::filesystem;

int pathDosentExists(const string_view &path) {
    cout << "Path : " << path << " n'existe pas ou n'est pas accessible.\n" ;
    return EXIT_FAILURE;
}

// --inputA <path> --output <path>
/*int program_option::parse_codon_count(const vector<string_view> &argv) {
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
    << "\t" << INPUTA << "\tChemin vers le dossier qui contient les fichiers contigs ou il faut compter les codons.\n\n"
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
}*/



