//
// Created by Florian Claisse on 02/05/2023.
//

#ifndef CONTIGDIFF_DIRECTORY_H
#define CONTIGDIFF_DIRECTORY_H

#include <iostream>
#include <sys/stat.h>

/**
 * Check if a file exists
 * @param name Path to the file
 * @return True if the file exists, false otherwise
 */
bool fileExists(const std::string &name);

std::string fileNameWithoutExtension(const std::string &filePath);
std::string fileName(const std::string &filePath);
std::string removeExtension(const std::string &value);

#endif //CONTIGDIFF_DIRECTORY_H
