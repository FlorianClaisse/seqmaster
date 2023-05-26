//
// Created by Florian Claisse on 22/05/2023.
//

#ifndef CONTIG_CENT_PERCENT_H
#define CONTIG_CENT_PERCENT_H

#include <map>
#include <string>

#include "contig_diff.h"
#include "parser.h"

namespace contig_diff::cent_percent {
    std::map<std::string, contig_diff::Common> main(const contig_diff::option &options);
}

#endif //CONTIG_CENT_PERCENT_H
