//
// Created by Florian Claisse on 02/05/2023.
//

#include <fstream>

#include "Headers/FastaDecoder.h"
#include "Headers/Directory.h"

void decodeFasta(std::string &inputPath) {
    if (!fileExists(inputPath)) {
        std::cerr << "File " << inputPath << " does not exist" << std::endl;
        exit(EXIT_FAILURE);
    }

    // open file
    std::ifstream inputFile;
    inputFile.open(inputPath);

    // open output file
    std::string outputPath(removeExtension(inputPath).append(".fastaline"));
    std::ofstream outputFile;
    outputFile.open(outputPath, std::ios::out | std::ios::trunc);

    // remove all /n inside all contig
    std::string lineRead;
    bool first(true);
    while (getline(inputFile, lineRead)) {
        if (lineRead.at(0) == '>') {
            if (!first) outputFile << std::endl;
            outputFile << lineRead << std::endl;
            if (first) first = false;
        } else {
            outputFile << lineRead;
        }
    }


    outputFile.close();
    inputFile.close();

}


