//
// Created by Florian Claisse on 02/05/2023.
//

#include <fstream>
#include <functional>
#ifdef __linux__
#include <algorithm>
#endif

#include "include/fasta.h"
#include "include/directory.h"

using namespace std;
namespace fs = std::filesystem;

int fasta::to_fastaline(const std::filesystem::path &filePath) {
    ifstream inputFile;
    inputFile.open(filePath);

    ofstream outputFile(fs::path(filePath).replace_extension("fastaline"), ios::trunc);

    string lineRead;
    bool first(true);
    while(getline(inputFile, lineRead)) {
        if (lineRead.at(0) == '>') {
            if (!first) outputFile << endl;
            outputFile << lineRead << endl;
            if (first) first = false;
        } else {
            transform(lineRead.begin(), lineRead.end(), lineRead.begin(), ::toupper);
            outputFile << lineRead;
        }
    }

    inputFile.close();
    outputFile.close();
    return EXIT_SUCCESS;
}

int fasta::directory_to_fasta_line(const std::filesystem::path &directoryPath) {
    for (const auto &currentFile : fs::directory_iterator(directoryPath)) {
        if (!directory::is_fasta_file(currentFile)) continue;
        if (directory::is_result_file(currentFile)) continue;
        to_fastaline(currentFile);
    }
    return EXIT_SUCCESS;
}

map<string, string> fasta::decode_fastaline(const fs::path &filePath) {
    ifstream inputFile(filePath);

    map<string, string> result;

    string contig, sequence;
    while(!inputFile.eof()) {
        getline(inputFile, contig);
        getline(inputFile, sequence);
        result[contig] = sequence;
    }
    inputFile.close();

    return result;
}

bool fasta::find_contig(const fs::path &filePath, const string &contig) {
    ifstream test_file(filePath);

    string lineRead;
    while(getline(test_file, lineRead)) {
        if (lineRead.at(0) == '>') continue;
        if (lineRead.find(contig) != string::npos) return true;
    }
    return false;
}

// 100%
void fasta::find_contig(const fs::path &file_path, const map<string, string> &contigs, bool nucleic, function<void(const string&, const string&, const string&)> func) {
    ifstream test_file(file_path);

    string line_read, name;
    while(getline(test_file, line_read)) {
        if (line_read.at(0) == '>') name = line_read;
        else {
            size_t pos = 0;
            for(const auto &contig : contigs) {
                while((pos = line_read.find(contig.second, pos)) != string::npos) {
                    if (nucleic) func(contig.first, name, contig.second);
                    else func(contig.first, name, line_read);
                    pos++;
                }
                pos = 0;
            }
        }
    }
}

// < 100%
void fasta::find_contigs(const fs::path &file_path, const map<string, string> &contigs, int maxErrorPercentage, bool nucleic, function<void(const string&, const string&, const string&, double)> func) {
    ifstream  test_file(file_path);

    string line_read, name;
    while(getline(test_file, line_read)) {
        if (line_read.at(0) == '>') name = line_read;
        else {
            unsigned long text_size(line_read.size());
            for (const auto &contig : contigs) {
                // equalSearch 
                unsigned long pattern_size(contig.second.size());
                unsigned long maxError((pattern_size * maxErrorPercentage / 100));
                unsigned long error;
                for (unsigned long i = 0; i < text_size; i++) {
                    error = 0;
                    if (i + pattern_size > text_size) {
                        error += ((pattern_size + i) - text_size);
                        if (error > maxError) break;
                    }
                    for (unsigned long j = 0; j < pattern_size; j++) {
                        if (line_read[i + j] != contig.second[j]) {
                            error++;
                            if (error > maxError) break;
                        }
                    }
                    if (error <= maxError) {
                        if (nucleic) func(contig.first, name, line_read.substr(i, pattern_size), (((double)error) / ((double)pattern_size)) * 100.0);
                        else func(contig.first, name, line_read, (((double)error) / ((double)pattern_size)) * 100.0);
                    }

                }   
            }
        }
    }
}


