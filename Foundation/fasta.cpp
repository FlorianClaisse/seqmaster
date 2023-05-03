//
// Created by Florian Claisse on 02/05/2023.
//

#include <fstream>

#include "Headers/fasta.h"
#include "Headers/directory.h"

using namespace std;
namespace fs = std::filesystem;

int fasta::to_fasta_line(const fs::path &filePath) {
    if (!fs::exists(filePath)) {
        cout << "Path : " << filePath << ", le fichier ou le dossier n'existe pas." << endl;
        return EXIT_FAILURE;
    }
    if (!fs::is_regular_file(filePath)) {
        cout << "Path : " << filePath << ", n'est pas un fichier." << endl;
        return EXIT_FAILURE;
    }

    ifstream inputFile;
    inputFile.open(filePath);

    fs::path outputPath(directory::removeExtension(filePath).append(".fastaline"));
    ofstream outputFile;
    outputFile.open(outputPath, ios::out | ios::trunc);

    string lineRead;
    bool first(true);
    while(getline(inputFile, lineRead)) {
        if (lineRead.at(0) == '>') {
            if (!first) outputFile << endl;
            outputFile << lineRead << endl;
            if (first) first = false;
        } else outputFile << lineRead;
    }

    inputFile.close();
    outputFile.close();
    return EXIT_SUCCESS;
}

bool fasta::find_contig(const fs::path &filePath, const string &contig) {
    ifstream test_file;
    test_file.open(filePath);

    string lineRead;
    while(getline(test_file, lineRead)) {
        if (lineRead.at(0) == '>') continue;
        if (lineRead.find(contig) != string::npos) return true;
    }
    return false;
}

map<string, bool> fasta::find_contigs(const fs::path &file_path,  const vector<string> &contigs) {
    ifstream test_file;
    test_file.open(file_path);

    map<string, bool> result;
    for (const auto &contig : contigs) {
        result[contig] = false;
    }

    string line_read;
    while(getline(test_file, line_read)) {
        if (line_read.at(0) == '>') continue;
        for(const auto &contig : contigs) {
            if (line_read.find(contig) != string::npos) {
                result[contig] = true;
            }
        }
    }

    return result;
}


