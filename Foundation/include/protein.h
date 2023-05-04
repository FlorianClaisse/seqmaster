#ifndef CONTIGDIFF_PROTEIN_H
#define CONTIGDIFF_PROTEIN_H

#include <string>
#include <map>

namespace protein {
    bool is_protein_file(const std::string &filePath);
    bool is_protein_line_file(const std::string &filePath);
    bool find_protein(const std::string &filePath, const std::string &protein);
    std::map<std::string, bool> find_proteins(const std::string &filePath, const std::vector<std::string> &proteins);
}

#endif