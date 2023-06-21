//
// Created by Florian Claisse on 02/05/2023.
//

#include <iostream>
#include <thread>

#include "include/find_all.hpp"
#include "include/search.hpp"

#include <seqan3/alignment/configuration/all.hpp>
#include <seqan3/io/sequence_file/all.hpp>

#include "../Utils/include/termcolor.hpp"
#include "../Utils/include/directory.h"
#include "../Utils/include/file.h"

#define NUCL "nucl"
#define PROT "prot"

using namespace std;
namespace fs = std::filesystem;

namespace find_all {
    int error_message(const string &message);
    int check_options(const fs::path &inputA, const fs::path &inputB, const fs::path &output, const string &type, int accept, int threads);
}

int find_all::error_message(const string &message) {
    cout << termcolor::red << termcolor::bold
         << message
         << termcolor::reset;
    return -1;
}

int find_all::check_options(const fs::path &inputA, const fs::path &inputB, const fs::path &output, const string &type, int accept, int threads) {
    if (threads <= 0) return error_message("You want to run the program on a number of threads less than or equal to 0, go try again ;).\n");
    if (type != PROT && type != NUCL) return error_message("For the type, the value must be \"prot\" or \"nucl\" did you read the docs ?.\n");
    if (accept < 0 || accept > 100) return error_message("You just gave an acceptance percentage lower than 0 or higher than 100. That's not very clever.\n");
    if (!file::is_fasta(inputA)) return error_message("You need to give fasta file in inputA. Did you really read the doc.\n");
    if (!is_directory(inputB)) return error_message("Input B is not a folder or does not exist. Try again the next one is the right one.\n");
    if(!directory::create_directories(output)) return error_message("Unable to find/create the output folder are you sure you have given a valid path.\n");

    for (const auto &filePath: fs::directory_iterator(output))
        fs::remove_all(filePath);

    return 0;
}

int find_all::main(const fs::path &inputA, const fs::path &inputB, const fs::path &output, const string &type, int accept, int threads) {

    if (check_options(inputA, inputB, output, type, accept, threads) != 0) return -1;

    file::clean(inputA);
    directory::clean_fastas(inputB);

    param options = {inputA, inputB, output, (type == NUCL), (100 - accept), threads};

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
        find_all::search<traits_t>(options, config);
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
        find_all::search<traits_t>(options, config);
    }
    return 0;
}