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

    template<typename sequence_t>
    std::vector<typename std::iterator_traits<sequence_t>::value_type>subsequence(sequence_t begin, sequence_t end) {
        using value_type = typename std::iterator_traits<sequence_t>::value_type;
        std::vector<value_type> subsequence;
        while(begin != end) {
            subsequence.push_back(*begin);
            ++begin;
        }
        return subsequence;
    }

    template<typename traits_t, typename config_t>
    void search_with_2_files(const std::filesystem::path &firstPath, const std::filesystem::path &secondPath, const contig_diff::param &options, const config_t &config) {
        seqan3::sequence_file_input<traits_t> f_in{firstPath}, s_in{secondPath};

        using record_t = decltype(f_in)::record_type;
        // Get all records for first and second file.
        std::vector<record_t> f_records, s_records;
        std::ranges::copy(f_in, std::back_inserter(f_records));
        std::ranges::copy(s_in, std::back_inserter(s_records));

        using sequence_t = decltype(f_in)::sequence_type;
        using pair_t = decltype(std::tie(f_records[0].sequence(), s_records[0].sequence()));
        // [ref_sequence, (id, sequence)] id = filename_contigname\tpercentage
        std::unordered_map<sequence_t, std::tuple<std::string, sequence_t>> common;
        for (const auto &record: f_records) {
            for (unsigned long i = record.sequence().size(); i >= 0; i++) {
                auto sub = subsequence(record.sequence().begin(), record.sequence().begin() + i);
                std::string test = sub;
                if (common.find(static_cast<const sequence_t&>(sub)) != common.end()) continue;

                std::vector<pair_t> source;
                for (long i = 0; i < s_records.size(); i++) {
                    source.push_back(std::tie(sub, s_records[i].sequence()));
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
                    double error = ((double)(100 * abs(best_score))) / s_records[index].sequence().size();
                    if (error <= options.error_rate) {
                        auto best_sub = subsequence(s_records[index].sequence().begin() + best_start, s_records[index].sequence().begin() + best_end);
                        //common[sub] = {(secondPath.filename().string() + "_" + s_records[index].id() + "\t" + std::to_string(100 - error)), best_sub};
                    }
                }
            }
        }

        //return common;
        //seqan3::debug_stream << common << '\n';
    }
}


#endif //CONTIG_CONTIG_DIFF_SEARCH_HPP
