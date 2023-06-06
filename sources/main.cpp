#include "Utils/include/parser.h"

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

using namespace seqan3::literals;

int main(int argc, const char* argv[]) {


    /*auto p1 = "ZZZZZZZZZZZZZZZZZZZZMKFTAALHRIIIGKKEVDETZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZ"_aa27;
    auto p2 = "MKCTAAVHRITIGKAEVDEH"_aa27;

    using traits_t = seqan3::sequence_file_input_default_traits_aa;
    using record_t = typename seqan3::sequence_file_input<traits_t>::record_type;
    using sequence_t = typename seqan3::sequence_file_input<traits_t>::sequence_type;
    using pair_t = typename std::tuple<sequence_t, sequence_t>;

    std::vector<pair_t> source;
    for (int i = 0; i < 20000000; i++)
        source.push_back(std::tie(p1, p2));

    seqan3::configuration const config = seqan3::align_cfg::method_global{
            seqan3::align_cfg::free_end_gaps_sequence1_leading{false},
            seqan3::align_cfg::free_end_gaps_sequence2_leading{true},
            seqan3::align_cfg::free_end_gaps_sequence1_trailing{false},
            seqan3::align_cfg::free_end_gaps_sequence2_trailing{true}}
                                         | seqan3::align_cfg::scoring_scheme{seqan3::aminoacid_scoring_scheme{
            seqan3::match_score{0},
            seqan3::mismatch_score{-1}}};

    auto res = seqan3::align_pairwise(source, config);*/
    //seqan3::debug_stream << res << '\n';
    parser::parse(argc, argv);

    return EXIT_SUCCESS;
}
