//
// Created by Florian Claisse on 02/05/2023.
//

#ifndef CONTIGDIFF_DIRECTORY_H
#define CONTIGDIFF_DIRECTORY_H

#include <iostream>
#include <vector>

namespace directory {
    bool have_extension(const std::string &filePath, const std::string &extension);
    bool have_extension(const std::string &filePath, const std::vector<std::string> &extensions);
    std::string fileNameWithoutExtension(const std::string &filePath);
    std::string fileName(const std::string &filePath);
    std::string removeExtension(const std::string &value);
}

#endif //CONTIGDIFF_DIRECTORY_H
