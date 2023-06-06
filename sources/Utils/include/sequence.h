//
// Created by Florian Claisse on 05/06/2023.
//

#ifndef CONTIG_SEQUENCE_H
#define CONTIG_SEQUENCE_H

#include <vector>

namespace sequence {
    template<typename sequence_t>
    std::vector<typename std::iterator_traits<sequence_t>::value_type>subsequence(sequence_t begin, sequence_t end) {
        return {begin, end};
    }
}

#endif //CONTIG_SEQUENCE_H
