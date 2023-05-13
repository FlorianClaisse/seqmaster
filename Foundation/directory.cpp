//
// Created by Florian Claisse on 02/05/2023.
//

#include "include/directory.h"

using namespace std;
namespace fs = std::filesystem;

bool foundation::have_extension(const fs::path &filePath, const vector<string> &extensions) {
    string file_extension = filePath.extension();
    for (const auto &extension: extensions) {
        if (file_extension == extension) return true;
    }

    return false;
}
