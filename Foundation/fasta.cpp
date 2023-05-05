//
// Created by Florian Claisse on 02/05/2023.
//

#include <fstream>
#include <tuple>
#include <functional>
#ifdef __linux__
#include <algorithm>
#endif

#include "include/fasta.h"
#include "include/directory.h"

using namespace std;
namespace fs = std::filesystem;

int fasta::to_fasta_line(const fs::path &filePath) {
    ifstream inputFile;
    inputFile.open(filePath);

    fs::path outputPath(directory::removeExtension(filePath).append(".fastaline"));
    ofstream outputFile;
    outputFile.open(outputPath, ios::trunc);

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

map<string, string> fasta::decode_fastaline(const fs::path &filePath) {
    ifstream inputFile;
    inputFile.open(filePath);

    map<string, string> result;

    string contig;
    string sequence;
    while(!inputFile.eof()) {
        getline(inputFile, contig);
        getline(inputFile, sequence);
        result[contig] = sequence;
    }
    inputFile.close();

    return result;
}

bool fasta::is_fasta_file(const fs::path &filePath) {
    if (!is_regular_file(filePath)) return false;

    vector<string> extensions = {"fasta", "fna", "faa", "ffn", "fa", "fas"};
    return directory::have_extension(filePath, extensions);
}

bool fasta::is_fastaline_file(const fs::path &filePath) {
    return is_regular_file(filePath) && directory::have_extension(filePath, "fastaline");
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

void fasta::find_contig(const fs::path &file_path, const map<string, string> &contigs, bool nucleic, function<void(const string&, const string&, const string&)> func) {
    ifstream test_file;
    test_file.open(file_path);

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

/*map<string, double> equalSearch(const string &text, const string &pattern, double maxErrorPercentage) {
    unsigned long text_size(text.size());
    unsigned long pattern_size(pattern.size());
    unsigned long maxError((unsigned long) (pattern_size * maxErrorPercentage / 100.0));
    unsigned long error(0);

    for (unsigned long i = 0; i < text_size; i++) {
        for (unsigned long j = 0; j < pattern_size; j++) {
            if (i + j > text_size) {
                error += (pattern_size - j);
                break;
            }
            if (text[i + j] != pattern[j]) error++;
            if (error >= result) break;
        }
        result = min(result, error);
        error = 0;
    }

    return (((double)result) / ((double)pattern_size)) * 100.0;
}*/

void fasta::find_contigs(const fs::path &file_path, const map<string, string> &contigs, int maxErrorPercentage, bool nucleic, function<void(const string&, const string&, const string&, double)> func) {
    ifstream  test_file;
    test_file.open(file_path);

    string line_read, name;
    while(getline(test_file, line_read)) {
        if (line_read.at(0) == '>') name = line_read;
        else {
            unsigned long text_size(line_read.size());
            for (const auto &contig : contigs) {
                // equalSearch 
                unsigned long pattern_size(contig.second.size());
                unsigned long maxError((unsigned long) (pattern_size * maxErrorPercentage / 100));
                unsigned long error(0);
                for (unsigned long i = 0; i < text_size; i++) {
                    for (unsigned long j = 0; j < pattern_size; j++) {
                        if (i + j > text_size) {
                            error += (pattern_size - j);
                            break;
                        }
                        if (line_read[i + j] != contig.second[j]) error++;
                        if (error > maxError) break;
                    }
                    if (error <= maxError) {
                        if (nucleic) func(contig.first, name, contig.second, (((double)error) / ((double)pattern_size)) * 100.0);
                        else func(contig.first, name, line_read, (((double)error) / ((double)pattern_size)) * 100.0);
                    }
                    error = 0;
                }   
            }
        }
    }
}


