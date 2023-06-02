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

#include <seqan3/io/sequence_file/all.hpp>
#include <seqan3/alignment/pairwise/align_pairwise.hpp>
#include <seqan3/alignment/configuration/all.hpp>
#include <seqan3/alphabet/all.hpp>
#include <seqan3/core/debug_stream.hpp>

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

    template<typename traits_t, typename config_t, typename out_t>
    void search_with_2_files(const path &firstPath, const path &secondPath, const param &options, const config_t &config, out_t &output) {
        seqan3::sequence_file_input<traits_t> f_in{firstPath}, s_in{secondPath};

        using record_t = decltype(f_in)::record_type;
        // Get all records for first and second file.
        vector<record_t> f_records, s_records;
        std::ranges::copy(f_in, std::back_inserter(f_records));
        std::ranges::copy(s_in, std::back_inserter(s_records));

        using sequence_t = decltype(f_in)::sequence_type;
        using pair_t = decltype(std::tie(f_records[0].sequence(), s_records[0].sequence()));
        // [(id, sequence)] id = filename_contigname\tpercentage
        vector<pair<string, sequence_t>> common;
        for (const auto &record: f_records) {
            for (unsigned long i = record.sequence().size(); i > 0; i--) {
                auto sub = subsequence(record.sequence().begin(), record.sequence().begin() + i);
                seqan3::debug_stream << sub << '\n';
                vector<pair_t> source;
                for (auto &s_record: s_records) {
                    source.push_back(std::tie(sub, s_record.sequence()));
                }
                long best_score{LONG_MIN}, index{-1}, best_start{-1}, best_end{-1};
                for (const auto &res : seqan3::align_pairwise(source, config)) {
                    if (res.score() > best_score) {
                        best_score = res.score();
                        index = res.sequence2_id();
                        best_start = res.sequence2_begin_position();
                        best_end = res.sequence2_end_position();
                    }
                }
                if (best_score != LONG_MIN) {
                    double error = ((double)(100 * abs(best_score))) / record.sequence().size();
                    if (error <= options.error_rate) {
                        auto best_sub = subsequence(s_records[index].sequence().begin() + best_start, s_records[index].sequence().begin() + best_end);
                        common.push_back({(secondPath.stem().string() + "_" + s_records[index].id() + "\t" + std::to_string(100 - error) + '%'), best_sub});
                        break;
                    }
                }
            }
        }

        for (const auto &value: common) {
            add_to(output, value.first, value.second);
        }
    }

    template<typename traits_t, typename config_t>
    void check_common(const path &input, const path &dir, const config_t &config) {
        seqan3::sequence_file_input<traits_t> f_in{input};

        using record_t = decltype(f_in)::record_type;
        vector<record_t> records;
        std::ranges::copy(f_in, std::back_inserter(records));

        for (const auto record: records) {
            for (const auto &test_path: std::filesystem::directory_iterator(dir)) {
                vector<record_t> test_records;
                std::ranges::copy(f_in, std::back_inserter(test_records));


            }
        }
    }
}


#endif //CONTIG_CONTIG_DIFF_SEARCH_HPP
