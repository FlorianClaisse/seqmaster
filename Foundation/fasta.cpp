//
// Created by Florian Claisse on 02/05/2023.
//

#include <fstream>

#include "Headers/fasta.h"
#include "Headers/directory.h"

int fasta::to_fasta_line(const std::filesystem::path &filePath) {
    if (!std::filesystem::exists(filePath)) {
        std::cout << "Path : " << filePath << ", le fichier ou le dossier n'existe pas." << std::endl;
        return EXIT_FAILURE;
    }
    if (!std::filesystem::is_regular_file(filePath)) {
        std::cout << "Path : " << filePath << ", n'est pas un fichier." << std::endl;
        return EXIT_FAILURE;
    }

    std::ifstream inputFile;
    inputFile.open(filePath);

    std::filesystem::path outputPath(directory::removeExtension(filePath).append(".fastaline"));
    std::ofstream outputFile;
    outputFile.open(outputPath, std::ios::out | std::ios::trunc);

    std::string lineRead;
    bool first(true);
    while(getline(inputFile, lineRead)) {
        if (lineRead.at(0) == '>') {
            if (!first) outputFile << std::endl;
            outputFile << lineRead << std::endl;
            if (first) first = false;
        } else outputFile << lineRead;
    }

    inputFile.close();
    outputFile.close();
    return EXIT_SUCCESS;
}

bool fasta::find_contig(const std::filesystem::path &filePath, const std::string &contig) {
    std::ifstream test_file;
    test_file.open(filePath);

    std::string lineRead;
    while(getline(test_file, lineRead)) {
        if (lineRead.at(0) == '>') continue;
        if (lineRead.find(contig) != std::string::npos) return true;
    }
    return false;
}

std::map<std::string, bool> fasta::find_contigs(const std::filesystem::path &file_path,  const std::vector<std::string> &contigs) {
    std::ifstream test_file;
    test_file.open(file_path);

    std::map<std::string, bool> result;
    for (const auto &contig : contigs) {
        result[contig] = false;
    }

    std::string line_read;
    while(getline(test_file, line_read)) {
        if (line_read.at(0) == '>') continue;
        for(const auto &contig : contigs) {
            if (line_read.find(contig) != std::string::npos) {
                result[contig] = true;
            }
        }
    }

    return result;
}


