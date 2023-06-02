#include "Utils/include/parser.h"

#include <filesystem>
#include <vector>
#include <ranges>

#include <seqan3/io/sequence_file/all.hpp>
#include <seqan3/alignment/pairwise/align_pairwise.hpp>
#include <seqan3/alignment/configuration/all.hpp>
#include <seqan3/alphabet/all.hpp>
#include <seqan3/core/debug_stream.hpp>

int main(int argc, const char* argv[]) {
    std::vector<seqan3::sequence_file_input<seqan3::sequence_file_input_default_traits_dna>::record_type> records;

    for (const auto &path: std::filesystem::directory_iterator(argv[1])) {
        seqan3::sequence_file_input<seqan3::sequence_file_input_default_traits_dna> f_in{path};
        std::ranges::copy(f_in, std::back_inserter(records));
    }

    //parser::parse(argc, argv);

    return EXIT_SUCCESS;
}
