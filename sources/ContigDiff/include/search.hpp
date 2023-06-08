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

// TODO: Ajouter les uniques d'un fichier
// TODO: optimiser pour de la recherche à 100%
// TODO: Ne pas chercher de duplicata
// TODO: Ne pas chercher de duplicata à l'intérieur d'un même fichier
// TODO: ...

namespace contig_diff {

    template<typename traits_t, typename config_t>
    class Search {
    private:
        using record_t = typename seqan3::sequence_file_input<traits_t>::record_type;
        using sequence_t = typename seqan3::sequence_file_input<traits_t>::sequence_type;
        using pair_t = typename std::tuple<sequence_t, sequence_t>;

    private:
        struct BestInfo {
            unsigned long fileIndex;
            unsigned long sequenceIndex;
            double error_percentage;
            unsigned long sequence_start;
            unsigned long sequence_end;
        };

    private:
        contig_diff::param options;
        config_t config;

    public:
        Search(contig_diff::param options, const config_t config): options(std::move(options)), config(std::move(config)) { }

    public:
        std::string search_common(const std::filesystem::path &dir) {
            directory::create_directories(options.output / ("Specific-" + dir.filename().string()));
            std::string commonPath{options.output / ("Common-" + dir.filename().string() + ".fasta")};
            seqan3::sequence_file_output common_out{commonPath};

            // [{filename, {record}}]
            std::vector<std::pair<std::string, std::vector<record_t>>> all_records;
            directory::decode_all_fasta<traits_t>(dir, all_records);

            std::set<sequence_t> global_common;
            for(unsigned long i = 0; i < all_records.size(); i++) { // Visit all file
                // [base_sequence, {file index, percentage}]
                std::map<sequence_t, BestInfo> current_common;
                for (unsigned long j = 0; j < all_records.size(); j++) { // Check inside all file
                    if (j == i) continue;

                    std::vector<std::vector<pair_t>> sources;
                    generate_test_source(all_records[i].second, current_common, all_records[j].second, sources, global_common);

                    for (auto source: sources) {
                        sequence_t sequence = (get<0>(source[0]));
                        long max_error = (sequence.size() * options.error_rate) / 100;

                        unsigned long best_nb_error{LONG_MAX}, best_index, best_start, best_end;
                        for (auto &res: seqan3::align_pairwise(source, config)) {
                            long nb_error = abs(res.score());

                            if (nb_error <= max_error && nb_error < best_nb_error) {
                                best_nb_error = nb_error;
                                best_index = res.sequence2_id();
                                best_start = res.sequence2_begin_position();
                                best_end = res.sequence2_end_position();
                            }
                        }

                        if (current_common.find(sequence) == current_common.end() && best_nb_error != LONG_MAX)
                            current_common[sequence] = {j, best_index, ((double)(100 * best_nb_error) / sequence.size()), best_start, best_end};
                        else if (best_nb_error != LONG_MAX) {
                            auto it = current_common.find(sequence);
                            auto error_percentage = ((double)(100 * best_nb_error)) / sequence.size();
                            if (error_percentage < it->second.error_percentage)
                                it->second = {j, best_index, error_percentage, best_start, best_end};
                        } else
                            current_common.erase(sequence);
                    }

                    if (current_common.empty()) break;
                }

                for (const auto &value: current_common) {
                    auto best_record = all_records[value.second.fileIndex];
                    std::string best_name = (best_record.first + " -> " + best_record.second[value.second.sequenceIndex].id() + "\t" + std::to_string(100 - value.second.error_percentage));
                    sequence_t best_sub = sequence::subsequence(best_record.second[value.second.sequenceIndex].sequence().begin() + value.second.sequence_start, best_record.second[value.second.sequenceIndex].sequence().begin() + value.second.sequence_end);
                    add_to(common_out, best_name, best_sub);
                    global_common.insert(value.first);
                }
            }

            return commonPath;
        }

        // TODO Chercher si les communs de A sont dans B
        void check_common(const std::filesystem::path &dir, const std::filesystem::path &commonPath) {
            seqan3::sequence_file_output ab_out{options.output / ("Common-" + options.inputA.filename().string() + "-" + options.inputB.filename().string() + ".fasta")};
            seqan3::sequence_file_output a_not_b_out{options.output / ("Common-" + options.inputA.filename().string() + "-not-" + options.inputB.filename().string() + ".fasta")};

            std::map<std::string, sequence_t> all_common;
            file::decode_fasta<traits_t>(commonPath, all_common);

            std::vector<std::pair<std::string, std::vector<record_t>>> all_records;
            directory::decode_all_fasta<traits_t>(dir, all_records);

            for (const auto &record: all_records) {

                std::vector<std::string> seq_to_remove;
                for (const auto &test_common: all_common) {
                    std::vector<pair_t> source;
                    generate_test_source(test_common.second, record.second, source);

                    long max_error = (test_common.second.size() * options.error_rate) / 100;
                    for (const auto &res: seqan3::align_pairwise(source, config)) {
                        long nb_error = abs(res.score());

                        if (nb_error <= max_error) {
                            seq_to_remove.push_back(test_common.first);
                            auto error_percentage = ((double)(100 * nb_error)) / test_common.second.size();
                            auto seq = record.second[res.sequence2_id()];
                            add_to(ab_out,
                                   (record.first + " -> " + seq.id() + "\t" + std::to_string(100 - error_percentage)),
                                   (sequence::subsequence(seq.sequence().begin() + res.sequence2_begin_position(), seq.sequence().begin() + res.sequence2_end_position()))
                                   );
                        }
                    }
                }

                for (const auto &remove: seq_to_remove)
                    all_common.erase(remove);
            }

            for (const auto common: all_common)
                add_to(a_not_b_out,common.first, common.second);
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

        bool already_found(const std::set<sequence_t> &common, const sequence_t &sub) const {
            auto it = common.find(sub);
            if (it != common.end()) return true;

            std::vector<pair_t> source;
            for (auto &value: common)
                source.push_back(std::tie(value, sub));

            long max_error = (sub.size() * options.error_rate) / 100;
            for (auto &res : seqan3::align_pairwise(source, config))
                if (abs(res.score()) <= max_error) return true;

            return false;
        }

        void generate_test_source(const std::vector<record_t> &records, const std::map<sequence_t, BestInfo> &current_common,
                                  const std::vector<record_t> &test_records, std::vector<std::vector<pair_t>> &sources,
                                  const std::set<sequence_t> &global_common) {
            if (current_common.empty()) {
                for (const auto &record: records) {
                    if (already_found(global_common, record.sequence())) continue;

                    std::vector<pair_t> source;
                    for (const auto &test_record: test_records)
                        source.push_back(std::tie(record.sequence(), test_record.sequence()));

                    sources.push_back(source);
                }
            } else {
                for (const auto &sequence: current_common) {
                    std::vector<pair_t> source;
                    for (const auto &test_record: test_records)
                        source.push_back(std::tie(sequence.first, test_record.sequence()));

                    sources.push_back(source);
                }
            }
        }

        inline void generate_test_source(const sequence_t &sequence, const std::vector<record_t> &records, std::vector<pair_t> &source) {
            for (const auto &record: records)
                source.push_back(std::tie(sequence, record.sequence()));
        }
    };
}


/*for (const auto record: all_records[i].second) {
                    if (already_found(record.sequence())) continue;

                    long max_error = (record.sequence().size() * options.error_rate) / 100;

                    double best_percentage{MAXFLOAT};
                    long best_start{-1}, best_end{-1};
                    std::pair<long, long> best_index;
                    for (unsigned long j = 0; j < all_records.size(); j++) {
                        if (j == i) continue;

                        std::vector<pair_t> source;
                        for (auto &test_record: all_records[j].second)
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
                        auto best_record = all_records[best_index.first].second[best_index.second];
                        auto best_sub = sequence::subsequence(best_record.sequence().begin() + best_start, best_record.sequence().begin() + best_end);
                        std::string best_name = (best_record.id() + '\t' + std::to_string(100 - best_percentage));
                        add_to(common_out, best_name, best_sub);
                        if (common.find(record.sequence().size()) == common.end()) common[record.sequence().size()] = {record.sequence()};
                        else common.at(record.sequence().size()).push_back(record.sequence());
                    }*/

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

#endif //CONTIG_CONTIG_DIFF_SEARCH_HPP
