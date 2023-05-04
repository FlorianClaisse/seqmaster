
#include "include/protein.h"

bool protein::is_protein_file(const std::string &filePath) {
    std::ifstream file(filePath);
    std::string line;
    while (std::getline(file, line)) {
        if (line[0] == '>') {
            return true;
        }
    }
    return false;
}