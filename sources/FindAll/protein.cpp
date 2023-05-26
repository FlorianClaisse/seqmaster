//
// Created by Florian Claisse on 26/05/2023.
//

#include <string>
#include <ranges>
#include <vector>
#include <unordered_map>

#include <seqan3/alignment/configuration/all.hpp>
#include <seqan3/alignment/pairwise/align_pairwise.hpp>
#include <seqan3/io/sequence_file/all.hpp>
#include <seqan3/core/debug_stream.hpp>

#include "include/protein.h"

#include "../Utils/include/file.h"

using namespace std;
namespace fs = std::filesystem;

using traits_type = seqan3::sequence_file_input_default_traits_aa;

using types = seqan3::type_list<vector<seqan3::aa27>, string, vector<seqan3::phred42>>;
using fields = seqan3::fields<seqan3::field::seq, seqan3::field::id, seqan3::field::qual>;
using record_type = seqan3::sequence_record<types, fields>;

using sequence_types = seqan3::type_list<vector<seqan3::aa27>, string>;
using sequence_fields = seqan3::fields<seqan3::field::seq, seqan3::field::id>;
using sequence_record_type = seqan3::sequence_record<sequence_types, sequence_fields>;

template<typename ConfigType>
unordered_map<string, double> search_in(const fs::path &test_path, const ConfigType &config, const vector<record_type> &records, double error_rate, const fs::path &output) {
    seqan3::sequence_file_output output1{output / test_path.filename()};
    seqan3::sequence_file_input<traits_type> test_in{test_path};

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
                results[record.id()] = (100 - error);
                sequence_record_type write_record{record.sequence(), test_records[index].id() + " -> " + record.id() + " -> " + to_string(100 - error) + "%"};
                output1.push_back(write_record);
            }
        }
    }
    seqan3::debug_stream << results << '\n';
    return results;
}

ofstream* start_output(const fs::path &output) {
    ofstream *fout = file::write_open(output, ios::trunc);

    (*fout) << "Filename\n";
    return fout;
}

void find_all::protein::search(const fs::path &inputA, const fs::path &inputB, const fs::path &output, int accept, int threads) {

    seqan3::sequence_file_input<traits_type> fin{inputA};

    vector<record_type> records;
    ranges::copy(fin, back_inserter(records));

    const int error_rate = 100 - accept;

    seqan3::configuration const config = seqan3::align_cfg::method_global{
                                            seqan3::align_cfg::free_end_gaps_sequence1_leading{false},
                                            seqan3::align_cfg::free_end_gaps_sequence2_leading{true},
                                            seqan3::align_cfg::free_end_gaps_sequence1_trailing{false},
                                            seqan3::align_cfg::free_end_gaps_sequence2_trailing{true}}
                                            | seqan3::align_cfg::scoring_scheme{seqan3::aminoacid_scoring_scheme{
                                                                                    seqan3::match_score{0},
                                                                                    seqan3::mismatch_score{-1}}};

    ofstream *fout = start_output(output / "output.txt");
    for (const auto &test_path: fs::directory_iterator(inputB)) {
        if (!file::is_fasta(test_path)) continue;

        unordered_map<string, double> results = search_in(test_path, config, records, error_rate, output);
        seqan3::debug_stream << results << '\n';
        if (!results.empty()) {
            (*fout) << test_path.path().filename();
            for (const auto &value: results) {
                (*fout) << "\t" << value.first << " -> " << value.second << "%";
            }
        }
    }

    file::write_close(fout);
}