//
// Created by Florian Claisse on 01/06/2023.
//

#ifndef CONTIG_CONTIG_DIFF_SEARCH_HPP
#define CONTIG_CONTIG_DIFF_SEARCH_HPP

#include <filesystem>
#include <vector>
#include <ranges>
#include <unordered_map>
#include <type_traits>
#include <limits.h>

#include <seqan3/io/sequence_file/all.hpp>
#include <seqan3/alignment/pairwise/align_pairwise.hpp>
#include <seqan3/alignment/configuration/all.hpp>
#include <seqan3/alphabet/all.hpp>
#include <seqan3/core/debug_stream.hpp>

#include "../../Utils/include/file.h"

#include "contig_diff.h"

namespace contig_diff {

    using path = std::filesystem::path;
    using string = std::string;
    using std::vector;
    using std::iterator_traits;
    using std::pair;

    template<typename sequence_t>
    vector<typename iterator_traits<sequence_t>::value_type>subsequence(sequence_t begin, sequence_t end) {
        return {begin, end};
    }

    template<typename sequence_t, typename out_t>
    void add_to(out_t &output, const string &contigName, const sequence_t &sequence) {
        using sequence_types = seqan3::type_list<sequence_t, string>;
        using sequence_fields = seqan3::fields<seqan3::field::seq, seqan3::field::id>;
        using sequence_record_type = seqan3::sequence_record<sequence_types, sequence_fields>;

        sequence_record_type output_record = {sequence, contigName};
        output.push_back(output_record);
    }

    template<typename traits_t, typename config_t>
    void search_common(const path &dir, const param &options, const config_t &config) {
        seqan3::sequence_file_output f_out{options.output / (dir.parent_path().filename().string() + "-common.fasta")};
        vector<vector<typename seqan3::sequence_file_input<traits_t>::record_type>> all_records;

        for (const auto &path: std::filesystem::directory_iterator(dir)) {
            if (!file::is_fasta(path)) continue;
            seqan3::sequence_file_input<traits_t> f_in{path};
            vector<typename seqan3::sequence_file_input<traits_t>::record_type> records;
            std::ranges::copy(f_in, std::back_inserter(records));
            all_records.push_back(records);
        }

        using pair_t = decltype(std::tie(all_records[0][0].sequence(), all_records[0][0].sequence()));
        for (unsigned long i = 0; i < all_records.size(); i++) { // Check all file
            for (const auto &record: all_records[i]) { // Check all sequence in each file
                for (unsigned long size = record.sequence().size(); size > 0; size--) { // check all size inside each sequence
                    auto sub = subsequence(record.sequence().begin(), record.sequence().begin() + size);
                    seqan3::debug_stream << sub << '\n';

                    long max_error = (sub.size() * options.error_rate) / 100;

                    float best_percentage{MAXFLOAT};
                    long best_start{-1}, best_end{-1};
                    pair<long, long> best_index;
                    for (unsigned long j = 0; j < all_records.size(); j++) { // Check subsequence in each files
                        if (j == i) continue;
                        vector<pair_t> source;
                        for (auto &test_record: all_records[j]) {
                            source.push_back(std::tie(sub, test_record.sequence()));
                        }

                        long best_nb_error{LONG_MAX}, current_best_index{-1}, current_best_start{-1}, current_best_end{-1};
                        for (const auto &res : seqan3::align_pairwise(source, config)) {
                            //seqan3::debug_stream << res << '\n';
                            long nb_error = abs(res.score());
                            if (nb_error <= max_error && nb_error < best_nb_error) {
                                seqan3::debug_stream << "Good job" << "\n";
                                best_nb_error = nb_error;
                                current_best_index = res.sequence2_id();
                                current_best_start = res.sequence2_begin_position();
                                current_best_end = res.sequence2_end_position();
                            }
                        }

                        if (best_nb_error != LONG_MAX) {

                            double error = ((double)(100 * best_nb_error)) / size;
                            if (error < best_percentage) {
                                seqan3::debug_stream << "i'm here" << "\n";
                                best_percentage = error;
                                best_start = current_best_start;
                                best_end = current_best_end;
                                best_index = {j, current_best_index};
                            }
                        } else {
                            seqan3::debug_stream << "Oh shit" << "\n";
                            best_percentage = -1;
                            break;
                        }
                    }

                    if (best_percentage != -1) { // Subsequence not found inside all file
                        seqan3::debug_stream << "Good job2 " << "\n";
                        auto best_record = all_records[best_index.first][best_index.second];
                        auto best_sub = subsequence(best_record.sequence().begin() + best_start, best_record.sequence().begin() + best_end);
                        string best_name = (best_record.id() + '\t' + std::to_string(100 - best_percentage));
                        add_to(f_out, best_name, best_sub);
                        break;
                    }
                    // else continue to search
                }
            }
        }
    }
}


#endif //CONTIG_CONTIG_DIFF_SEARCH_HPP
