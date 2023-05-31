//
// Created by Florian Claisse on 31/05/2023.
//

#ifndef CONTIG_FINDALL_SEARCH_HPP
#define CONTIG_FINDALL_SEARCH_HPP


#include <string>
#include <filesystem>
#include <vector>
#include <unordered_map>
#include <ranges>

#include <seqan3/io/sequence_file/all.hpp>
#include <seqan3/alignment/pairwise/align_pairwise.hpp>
#include <seqan3/alphabet/all.hpp>

#include "../../Utils/include/file.h"

#include "find_all.hpp"

namespace find_all {
    template<typename sequence_t>
    std::vector<typename std::iterator_traits<sequence_t>::value_type>get_subsequence(sequence_t begin, sequence_t end) {
        using value_type = typename std::iterator_traits<sequence_t>::value_type;
        std::vector<value_type> subsequence;
        while(begin != end) {
            subsequence.push_back(*begin);
            ++begin;
        }
        return subsequence;
    }

    template<typename out_t, typename sequence_t>
    void generate_output(out_t &output, const std::string &contigName, const sequence_t &sequence) {
        using sequence_types = seqan3::type_list<sequence_t, std::string>;
        using sequence_fields = seqan3::fields<seqan3::field::seq, seqan3::field::id>;
        using sequence_record_type = seqan3::sequence_record<sequence_types, sequence_fields>;

        sequence_record_type output_record = {sequence, contigName};
        output.push_back(output_record);
    }

    template<typename traits_t, typename config_t, typename record_t>
    std::unordered_map<std::string, double> search_in(const std::filesystem::path &path, const std::vector<record_t> &records, const param &options, const config_t &config) {
        seqan3::sequence_file_output output{options.output / path.filename()};
        seqan3::sequence_file_input<traits_t> test_in{path};

        std::vector<record_t> test_records;
        std::ranges::copy(test_in, back_inserter(test_records));

        std::unordered_map<std::string, double> results;
        for (const auto &record : records) {
            long best_score{LONG_MIN}, index{-1}, best_start{-1}, best_end{-1};
            for (long i = 0; i < test_records.size(); i++) {
                auto result = seqan3::align_pairwise(tie(record.sequence(), test_records[i].sequence()), config);
                auto &res = *result.begin();
                if (res.score() > best_score) {
                    best_score = res.score();
                    index = i;
                    best_start = res.sequence2_begin_position();
                    best_end = res.sequence2_end_position();
                    if (best_score == 0) break;
                }
            }
            if (best_score != LONG_MIN) {
                double error = ((double)(100 * abs(best_score))) / test_records[index].sequence().size();
                if (error <= options.error_rate) {
                    results[record.id()] = (100 - error);
                    std::string contigName{test_records[index].id() + " -> " + record.id() + " -> " + std::to_string(100 - error) + "%"};
                    if (options.nucl) generate_output(output, contigName, get_subsequence(test_records[index].sequence().begin() + best_start, test_records[index].sequence().begin() + best_end));
                    else generate_output(output, contigName, test_records[index].sequence());
                }
            }
        }

        return results;
    }

    std::ofstream* start_output(const std::filesystem::path &outputPath) {
        std::ofstream *fout = file::write_open(outputPath, std::ios::trunc);

        (*fout) << "Filename\n";
        return fout;
    }

    template<typename traits_t, typename config_t>
    int search(const find_all::param &options, const config_t &config) {
        seqan3::sequence_file_input<traits_t> fin{options.inputA};

        using record_type = decltype(fin)::record_type;
        std::vector<record_type> records;
        std::ranges::copy(fin, back_inserter(records));

        std::ofstream *fout = start_output(options.output / "output.txt");
        for (const auto &test_path : std::filesystem::directory_iterator(options.inputB)) {
            if (!file::is_fasta(test_path)) continue;

            std::unordered_map<std::string, double> results = search_in<traits_t>(test_path, records, options, config);
            if (!results.empty()) {
                (*fout) << test_path.path().filename();
                for (const auto &value : results) {
                    (*fout) << "\t" << value.first << " -> " << value.second << "%";
                }
                (*fout) << "\n";
            }
        }

        file::write_close(fout);
        return 0;
    }
}

#endif //CONTIG_FINDALL_SEARCH_HPP
