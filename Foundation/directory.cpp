//
// Created by Florian Claisse on 02/05/2023.
//

#include "Headers/directory.h"

using namespace std;

bool directory::have_extension(const string &filePath, const string &extension) {
    return filePath.substr(filePath.find_last_of('.') + 1) == extension;
}

string directory::fileNameWithoutExtension(const string &filePath) {
    return removeExtension(fileName(filePath));
}

string directory::fileName(const string &filePath) {
    return filePath.substr(filePath.find_last_of('/') + 1);
}

string directory::removeExtension(const string &value) {
    string::size_type const p(value.find_last_of('.'));
    return value.substr(0, p);
}
