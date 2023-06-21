//
// Created by Florian Claisse on 10/05/2023.
//

#include <cstdint>
#include <iostream>
#include <filesystem>

#include "include/contig_diff.h"
#include "include/search.hpp"

#include "../Utils/include/termcolor.hpp"
#include "../Utils/include/directory.h"
#include "../Utils/include/file.h"

#define PROT "prot"
#define NUCL "nucl"

using namespace std;
namespace fs = std::filesystem;

namespace contig_diff {
    int error_message(const string &message);
    int check_options(const fs::path &inputA, const fs::path &inputB, const fs::path &output, const string &type, int accept, int threads);
}

int contig_diff::error_message(const string &message) {
    cout << termcolor::red << termcolor::bold
         << message
         << termcolor::reset;
    return -1;
}

int contig_diff::check_options(const fs::path &inputA, const fs::path &inputB, const fs::path &output, const string &type, int accept, int threads) {
    if (threads <= 0) return error_message("You want to run the program on a number of threads less than or equal to 0, go try again ;).\n");
    if (type != PROT && type != NUCL) return error_message("For the type, the value must be \"prot\" or \"nucl\" did you read the docs ?.\n");
    if (accept < 0 || accept > 100) return error_message("You just gave an acceptance percentage lower than 0 or higher than 100. That's not very clever.\n");
    if (!is_directory(inputA)) return error_message("Input A is not a folder or does not exist. Try again the next one is the right one.\n");
    if (!is_directory(inputB)) return error_message("Input B is not a folder or does not exist. Try again the next one is the right one.\n");
    if (!directory::create_directories(output)) return error_message("Unable to find/create the output folder are you sure you have given a valid path.\n");

    return 0;
}

template<typename traits_t, typename config_t>
int start_search(const contig_diff::param &options, const config_t &config) {
    
    contig_diff::Search<traits_t, config_t> search(options, config);

    cout << "Start all common in A.\n";
    std::string commonPath = search.search_common(true);

    cout << "\nCheck all common inside B\n";
    search.check_common(true, commonPath);

    cout << "\nStart all common inside B\n";
    search.search_common(false);

    cout << "\nFind unique inside all A files\n";
    std::string AuniquePath = search.search_unique(true);

    cout << "\nFind specific of A\n";
    search.search_specific(AuniquePath, false);

    cout << "\nFind unique inside all B files\n";
    std::string BuniquePath = search.search_unique(false);

    cout << "\nFind specific of B\n";
    search.search_specific(BuniquePath, true);

    return 0;
}

int contig_diff::main(const fs::path &inputA, const fs::path &inputB, const fs::path &output, const string &type, int accept, int threads, unsigned long min_size) {

    fs::remove_all(output);

    if (check_options(inputA, inputB, output, type, accept, threads) != 0) return -1;

    directory::clean_fastas(inputA);
    directory::clean_fastas(inputB);

    int nb_fasta = directory::count_fasta_file(inputA);
    if (nb_fasta < 2) {
        cout << "Il faut au minimum deux fichiers de type fasta dans le dossier A.\n";
        return -1;
    }
    nb_fasta = directory::count_fasta_file(inputB);
    if (nb_fasta < 1) {
        cout << "Il faut au minimum 1 fichier de type fasta dans le dossier B.\n";
        return -1;
    }

    param options = {inputA, inputB, output, (type == NUCL), (100 - accept), threads, min_size};

    if (options.nucl) {
        seqan3::configuration const config = seqan3::align_cfg::method_global{
                seqan3::align_cfg::free_end_gaps_sequence1_leading{false},
                seqan3::align_cfg::free_end_gaps_sequence2_leading{true},
                seqan3::align_cfg::free_end_gaps_sequence1_trailing{false},
                seqan3::align_cfg::free_end_gaps_sequence2_trailing{true}}
                                             | seqan3::align_cfg::scoring_scheme{seqan3::nucleotide_scoring_scheme{
                seqan3::match_score{0},
                seqan3::mismatch_score{-1}}}
        | seqan3::align_cfg::parallel{static_cast<uint32_t>(options.threads)};

        using traits_t = seqan3::sequence_file_input_default_traits_dna;

        start_search<traits_t>(options, config);
    } else {
        seqan3::configuration const config = seqan3::align_cfg::method_global{
                seqan3::align_cfg::free_end_gaps_sequence1_leading{false},
                seqan3::align_cfg::free_end_gaps_sequence2_leading{true},
                seqan3::align_cfg::free_end_gaps_sequence1_trailing{false},
                seqan3::align_cfg::free_end_gaps_sequence2_trailing{true}}
                                             | seqan3::align_cfg::scoring_scheme{seqan3::aminoacid_scoring_scheme{
                seqan3::match_score{0},
                seqan3::mismatch_score{-1}}}
        | seqan3::align_cfg::parallel{static_cast<uint32_t>(options.threads)};

        using traits_t = seqan3::sequence_file_input_default_traits_aa;

        start_search<traits_t>(options, config);
    }

    return 0;
}



