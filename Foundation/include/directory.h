//
// Created by Florian Claisse on 02/05/2023.
//

#ifndef CONTIGDIFF_DIRECTORY_H
#define CONTIGDIFF_DIRECTORY_H

#include <string>
#include <filesystem>
#include <vector>

namespace foundation {
    bool have_extension(const std::filesystem::path &filePath, const std::vector<std::string> &extensions);
}

#endif //CONTIGDIFF_DIRECTORY_H
