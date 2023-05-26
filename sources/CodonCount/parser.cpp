#include <iostream>

#include "include/parser.h"
#include "include/condo_count.h"

#define INPUT "--input"
#define OUTPUT "--output"

using namespace std;

namespace codon_count {
    int check_option(const codon_count::option &options);
    int usage();
}

// --input <path> --output <path>
int codon_count::parse(const std::vector <std::string_view> &argv) {
    if (argv.size() == 4) {
        if (argv[0] == INPUT && argv[2] == OUTPUT) {
            codon_count::option options {argv[1], argv[3]};
            if (check_option(options) == EXIT_FAILURE) return EXIT_FAILURE;

            //return codon_count::main(options);
        } else return usage();
    } else return usage();
}

int codon_count::check_option(const codon_count::option &options) {
    if (!is_directory(options.input)) {
        cout << "Path : " << options.input << " n'est pas un dossier.\n";
        return EXIT_FAILURE;
    }
    if (!is_directory(options.output)) {
        cout << "Path : " << options.output << " n'est pas un dossier.\n";
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}

int codon_count::usage() {
    cout << "Codon Count\n"
         << "Usage :\n"
         << "./Contig " << "--codoncount" << " " << INPUT << " <path> " << OUTPUT << " <path>\n\n"
         << "\t" << INPUT << "\t\tChemin vers le dossier qui contient les fichiers contigs ou il faut compter les codons.\n\n"
         << "\t" << OUTPUT << "\tChemin vers le dossier de sortie.\n";

    return EXIT_SUCCESS;
}