//
// Created by Florian Claisse on 01/06/2023.
//

#ifndef CONTIG_CONTIG_DIFF_SEARCH_HPP
#define CONTIG_CONTIG_DIFF_SEARCH_HPP

#include <filesystem>
#include <utility>
#include <vector>
#include <ranges>
#include <unordered_map>
#include <type_traits>

#include <seqan3/io/sequence_file/all.hpp>
#include <seqan3/alignment/pairwise/align_pairwise.hpp>
#include <seqan3/alignment/configuration/all.hpp>
#include <seqan3/alphabet/all.hpp>
#include <seqan3/core/debug_stream.hpp>

#include "../../Utils/include/file.h"
#include "../../Utils/include/directory.h"
#include "../../Utils/include/sequence.h"

#include "contig_diff.h"

namespace contig_diff {

    template<typename traits_t, typename config_t>
    class Search {
    private:
        using record_t = typename seqan3::sequence_file_input<traits_t>::record_type;
        using sequence_t = typename seqan3::sequence_file_input<traits_t>::sequence_type;
        using pair_t = typename std::tuple<sequence_t, sequence_t>;

    private:
        contig_diff::param options;
        config_t config;
        // [sequence_size, [original_sequence]]
        std::unordered_map<unsigned long, std::vector<sequence_t>> common;

    public:
        Search(contig_diff::param options, const config_t config): options(std::move(options)), config(std::move(config)), common() { }

    public:
        void search_common(const std::filesystem::path &dir) {
            directory::create_directories(options.output / (dir.filename()));

            seqan3::sequence_file_output common_out{options.output / (dir.filename().string() + "-common.fasta")};

            std::vector<std::vector<record_t>> all_records;
            directory::decode_all_fasta<traits_t>(dir, all_records);

            for(unsigned long i = 0; i < all_records.size(); i++) {
                long contig_counter{0};
                for (const auto record: all_records[i]) {
                    show_progress((i+1), all_records.size(), ++contig_counter, all_records[i].size());

                    long max_error = (record.sequence().size() * options.error_rate) / 100;

                    double best_percentage{MAXFLOAT};
                    long best_start{-1}, best_end{-1};
                    std::pair<long, long> best_index;
                    for (unsigned long j = 0; j < all_records.size(); j++) {
                        if (j == i) continue;

                        std::vector<pair_t> source;
                        for (auto &test_record: all_records[j])
                            source.push_back(std::tie(record.sequence(), test_record.sequence()));

                        long best_nb_error{LONG_MAX}, current_best_index{-1}, current_best_start{-1}, current_best_end{-1};
                        for (const auto &res: seqan3::align_pairwise(source, config)) {
                            long nb_error = abs(res.score());

                            if (nb_error <= max_error && nb_error < best_nb_error) {
                                best_nb_error = nb_error;
                                current_best_index = res.sequence2_id();
                                current_best_start = res.sequence2_begin_position();
                                current_best_end = res.sequence2_end_position();
                            }
                        }

                        if (best_nb_error != LONG_MAX) {
                            double error = ((double)(100 * best_nb_error)) / record.sequence().size();
                            if (error < best_percentage) {
                                best_percentage = error;
                                best_start = current_best_start;
                                best_end = current_best_end;
                                best_index = {j, current_best_index};
                            }
                        } else {
                            best_percentage = -1;
                            break;
                        }
                    }

                    if (best_percentage != -1) { // Subsequence found inside all file
                        auto best_record = all_records[best_index.first][best_index.second];
                        auto best_sub = sequence::subsequence(best_record.sequence().begin() + best_start, best_record.sequence().begin() + best_end);
                        std::string best_name = (best_record.id() + '\t' + std::to_string(100 - best_percentage));
                        add_to(common_out, best_name, best_sub);
                        if (common.find(record.sequence().size()) == common.end()) common[record.sequence().size()] = {record.sequence()};
                        else common.at(record.sequence().size()).push_back(record.sequence());
                    }
                }
            }
        }

        /*void search_common(const std::filesystem::path &dir) {
            seqan3::sequence_file_output f_out{options.output / (dir.filename().string() + "-common.fasta")};
            std::vector<std::vector<typename seqan3::sequence_file_input<traits_t>::record_type>> all_records;

            directory::decode_all_fasta<traits_t>(dir, all_records);

            for (unsigned long i = 0; i < all_records.size(); i++) { // Check all file
                long contig_counter{0};
                for (const auto &record: all_records[i]) { // Check all sequence in each file
                    show_progress((i+1), all_records.size(), ++contig_counter, all_records[i].size());

                    for (unsigned long size = record.sequence().size(); size > 0; size--) { // check all size inside each sequence
                        auto sub = sequence::subsequence(record.sequence().begin(), record.sequence().begin() + size);
                        if (already_exists(sub)) continue;

                        long max_error = (sub.size() * options.error_rate) / 100;

                        double best_percentage{MAXFLOAT};
                        long best_start{-1}, best_end{-1};
                        std::pair<long, long> best_index;
                        for (unsigned long j = 0; j < all_records.size(); j++) { // Check subsequence in each files
                            if (j == i) continue;
                            std::vector<pair_t> source;
                            for (auto &test_record: all_records[j])
                                source.push_back(std::tie(sub, test_record.sequence()));

                            long best_nb_error{LONG_MAX}, current_best_index{-1}, current_best_start{-1}, current_best_end{-1};
                            for (const auto &res : seqan3::align_pairwise(source, config)) {
                                //seqan3::debug_stream << res << '\n';
                                long nb_error = abs(res.score());
                                if (nb_error <= max_error && nb_error < best_nb_error) {
                                    best_nb_error = nb_error;
                                    current_best_index = res.sequence2_id();
                                    current_best_start = res.sequence2_begin_position();
                                    current_best_end = res.sequence2_end_position();
                                }
                            }

                            if (best_nb_error != LONG_MAX) {
                                double error = ((double)(100 * best_nb_error)) / size;
                                if (error < best_percentage) {
                                    best_percentage = error;
                                    best_start = current_best_start;
                                    best_end = current_best_end;
                                    best_index = {j, current_best_index};
                                }
                            } else {
                                best_percentage = -1;
                                break;
                            }
                        }

                        if (best_percentage != -1) { // Subsequence found inside all file
                            auto best_record = all_records[best_index.first][best_index.second];
                            auto best_sub = sequence::subsequence(best_record.sequence().begin() + best_start, best_record.sequence().begin() + best_end);
                            std::string best_name = (best_record.id() + '\t' + std::to_string(100 - best_percentage));
                            add_to(f_out, best_name, best_sub);
                            if (common.find(sub.size()) == common.end()) common[sub.size()] = {sub};
                            else common.at(sub.size()).push_back(sub);
                            break;
                        }
                        // else continue to search
                    }
                }
            }
        }*/

        // TODO Chercher si les communs de A sont dans B
        void check_common(const std::filesystem::path &dir) {

        }

    private:
        void show_progress(int currentFile, int totalFile, int currentRecord, int totalRecord) const {
            std::cout << "\r\033[K" << "File : " << currentFile << "/" << totalFile << " Finish at : " << (currentRecord * 100 / totalRecord) << "%";
            std::cout.flush();
        }

        template<typename out_t>
        void add_to(out_t &output, const std::string &contigName, const sequence_t &sequence) const {
            using sequence_types = seqan3::type_list<sequence_t, std::string>;
            using sequence_fields = seqan3::fields<seqan3::field::seq, seqan3::field::id>;
            using sequence_record_type = seqan3::sequence_record<sequence_types, sequence_fields>;

            sequence_record_type output_record = {sequence, contigName};
            output.push_back(output_record);
        }

        bool already_exists(const sequence_t &sub) const {
            auto it = common.find(sub.size());
            if (it == common.end()) return false;

            std::vector<pair_t> source;
            for (const auto &value: it->second)
                source.push_back(std::tie(value, sub));

            long max_error = (sub.size() * options.error_rate) / 100;
            for (const auto &res : seqan3::align_pairwise(source, config))
                if (abs(res.score()) <= max_error) return true;

            return false;
        }
    };
    // TODO Pour 100% pas besoin de faire un alignement juste de faire seq1 == seq2
}


#endif //CONTIG_CONTIG_DIFF_SEARCH_HPP
