//
// Created by Florian Claisse on 02/05/2023.
//

#include "Headers/directory.h"

#include <sys/stat.h>

bool directory::have_extension(const std::string &filePath, const std::string &extension) {
    return filePath.substr(filePath.find_last_of('.') + 1) == extension;
}

std::string directory::fileNameWithoutExtension(const std::string &filePath) {
    return removeExtension(fileName(filePath));
}

std::string directory::fileName(const std::string &filePath) {
    return filePath.substr(filePath.find_last_of('/') + 1);
}

std::string directory::removeExtension(const std::string &value) {
    std::string::size_type const p(value.find_last_of('.'));
    return value.substr(0, p);
}
