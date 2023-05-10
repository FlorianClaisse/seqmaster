//
// Created by Florian Claisse on 10/05/2023.
//

#ifndef CONTIG_CONTIG_DIFF_H
#define CONTIG_CONTIG_DIFF_H

#include "../Foundation/include/program_option.h"

namespace contig_diff {
    int start(const program_option::ContigDiff &options);

    typedef struct {
        unsigned long start_pos;
        std::string value;
    } Word;
}

#endif //CONTIG_CONTIG_DIFF_H
