//
// Created by Florian Claisse on 25/05/2023.
//

#ifndef CONTIG_FILE_H
#define CONTIG_FILE_H

#include <string>
#include <fstream>
#include <filesystem>
#include <vector>

#include <seqan3/io/sequence_file/all.hpp>


namespace file {
    std::ofstream* write_open(const std::filesystem::path &path, std::ios_base::openmode mode);
    std::ifstream* read_open(const std::filesystem::path &path);

    void write_close(std::ofstream *path);
    void read_close(std::ifstream *path);

    bool have_extension(const std::filesystem::path &path, const std::string &ext);
    bool have_extension(const std::filesystem::path &path, const std::vector<std::string> &exts);
    bool is_fasta(const std::filesystem::path &path);

    int to_fastaline(const std::filesystem::path &filePath);

    template<typename traits_t, typename record_t>
    void decode_fasta(const std::filesystem::path &path, std::vector<record_t> &all_records) {
        seqan3::sequence_file_input<traits_t> f_in{path};
        std::ranges::copy(f_in, std::back_inserter(all_records));
    }

    template<typename traits_t, typename sequence_t>
    void decode_fasta(const std::filesystem::path &path, std::map<std::string, sequence_t> &all_records) {
        seqan3::sequence_file_input<traits_t> f_in{path};

        for (auto record: f_in)
            all_records[record.id()] = record.sequence();
    }

    template<typename out_t, typename sequence_t>
    void add_to(out_t &output, const std::string &contigName, const sequence_t &sequence) {
        using sequence_types = seqan3::type_list<sequence_t, std::string>;
        using sequence_fields = seqan3::fields<seqan3::field::seq, seqan3::field::id>;
        using sequence_record_type = seqan3::sequence_record<sequence_types, sequence_fields>;

        sequence_record_type output_record = {sequence, contigName};
        output.push_back(output_record);
    }
}

#endif //CONTIG_FILE_H
