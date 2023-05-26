//
// Created by Florian Claisse on 02/05/2023.
//

#include <iostream>
#include <filesystem>
#include <ranges>
#include <vector>

#include <seqan3/alphabet/all.hpp>
#include <seqan3/alignment/configuration/all.hpp>
#include <seqan3/alignment/pairwise/align_pairwise.hpp>
#include <seqan3/io/sequence_file/all.hpp>
#include <seqan3/core/debug_stream.hpp>

#include "include/find_all.h"
#include "include/protein.h"

#include "../Utils/include/termcolor.hpp"
#include "../Utils/include/file.h"
#include "../Utils/include/directory.h"

#define NUCL "nucl"
#define PROT "prot"

using namespace std;
namespace fs = std::filesystem;

auto const convert_to_char_view = views::transform(
        [](auto const alph) {
            return seqan3::to_char(alph);
        });

int error_message(const string &message) {
    cout << termcolor::red << termcolor::bold
         << message
         << termcolor::reset;
    return -1;
}

int check_options(const fs::path &inputA, const fs::path &inputB, const fs::path &output, const string &type, int accept,
                  int threads) {
    if (threads <= 0) return error_message("You want to run the program on a number of threads less than or equal to 0, go try again ;).\n");
    if (type != PROT && type != NUCL) return error_message("For the type, the value must be \"prit\" or \"nucl\" did you read the docs ?.\n");
    if (accept <= 0 || accept >= 100) return error_message("You just gave an acceptance percentage lower than 0 or higher than 100. That's not very clever.\n");
    if (!file::is_fasta(inputA)) return error_message("You need to give fasta file in inputA. Did you really read the doc.\n");
    if (!is_directory(inputB)) return error_message("Input B is not a folder or does not exist. Try again the next one is the right one.\n");
    if(!directory::create_directories(output)) return error_message("Unable to find/create the output folder are you sure you have given a valid path.\n");

    return 0;
}

/*template<typename TraitsType, typename RecordType, typename Config>
void search_algo(const vector<RecordType> &records, const Config &config, const fs::path &test_path, double error_rate, const fs::path &output_path) {
    ofstream *output{directory::write_open(output_path, ios::trunc)};
    seqan3::sequence_file_input<TraitsType> test_in{test_path};

    using record_type = decltype(test_in)::record_type;
    vector<record_type> test_records;
    ranges::copy(test_in, back_inserter(test_records));

    unordered_map<string, double> results;
    for (const auto &record: records) {
        long best_score{LONG_MIN}, index{-1};
        for (long i = 0; i < test_records.size(); i++) {
            auto result = seqan3::align_pairwise(std::tie(record.sequence(), test_records[i].sequence()), config);
            auto &res = *result.begin();
            if (res.score() > best_score) {
                best_score = res.score();
                index = i;
            }
        }
        if (best_score != LONG_MIN) {
            double error = (double)(100 * abs(best_score)) / test_records[index].sequence().size();
            if (error <= error_rate) {
                results[record.id()] = error;
                (*output) << test_records[index].id() << " -> " << record.id() << " -> " << error << endl;
                auto v = record.sequence() | convert_to_char_view;
                (*output) << v << endl;
            }
        }
    }
}*/

int find_all::main(const fs::path &inputA, const fs::path &inputB, const fs::path &output, const string &type,
                   int accept, int threads) {

    check_options(inputA, inputB, output, type, accept, threads);

    protein::search(inputA, inputB, output, accept, threads);
    return 0;
}