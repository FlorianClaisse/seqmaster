#include "Utils/include/parser.h"
#include "Utils/include/termcolor.hpp"

#include <utility>

#include <seqan3/alignment/pairwise/align_pairwise.hpp>
#include <seqan3/alignment/scoring/all.hpp>
#include <seqan3/alphabet/nucleotide/all.hpp>
#include <seqan3/core/debug_stream.hpp>

using namespace std;

int main(int argc, const char* argv[]) {

    using namespace seqan3::literals;

    seqan3::dna4_vector s1 = "ATGTTCAG"_dna4;
    seqan3::dna4_vector s2 = "CTGAACAT"_dna4;

    // Configure the alignment kernel.
    auto config = seqan3::align_cfg::method_global{}
                  | seqan3::align_cfg::scoring_scheme{
            seqan3::nucleotide_scoring_scheme{seqan3::match_score{0}, seqan3::mismatch_score{-1}}};


    // Invoke the pairwise alignment which returns a lazy range over alignment results.
    auto results = seqan3::align_pairwise(std::tie(s1, s2), config);
    auto & res = *results.begin();
    seqan3::debug_stream << "Score: " << res << '\n';

    parser::parse(argc, argv);

    return EXIT_SUCCESS;
}
