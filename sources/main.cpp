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
    parser::parse(argc, argv);

    return EXIT_SUCCESS;
}
