//
// Created by Florian Claisse on 02/05/2023.
//

#include "Headers/Directory.h"

#include <sys/stat.h>

bool fileExists(const std::string &name) {
    struct stat buffer;
    return (stat(name.c_str(), &buffer) == 0);
}

std::string fileNameWithoutExtension(const std::string &filePath) {
    return removeExtension(fileName(filePath));
}

std::string fileName(const std::string &filePath) {
    return filePath.substr(filePath.find_last_of('/') + 1);
}

std::string removeExtension(const std::string &value) {
    std::string::size_type const p(value.find_last_of('.'));
    return value.substr(0, p);
}
