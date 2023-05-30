#include "Utils/include/parser.h"

#include <seqan3/alignment/pairwise/align_pairwise.hpp>
#include <seqan3/alignment/scoring/nucleotide_scoring_scheme.hpp>
#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/core/debug_stream.hpp>

int main(int argc, const char* argv[]) {
    
    parser::parse(argc, argv);

    return EXIT_SUCCESS;
}
