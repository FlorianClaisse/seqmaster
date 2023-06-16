//
// Created by Florian Claisse on 04/05/2023.
//

#include <map>
#include <iostream>
#include <fstream>
#include <algorithm>

#include "include/condo_count.h"

#include "../Utils/include/directory.h"
#include "../Utils/include/protein.h"

using namespace std;
namespace fs = std::filesystem;

namespace codon_count {
    void init_codon(map<string, int> &codon);
    void reset_codon(map<string, int> &codon);
    void init_output(ofstream *output);
    void add_to(ofstream *output, const map<string, int> &codon, const string &protName, int total, long gcNum, double gcpercent);
    int check_args(const fs::path &input, const fs::path &output);
}

void codon_count::init_output(std::ofstream *output) {
    for (const auto &value: protein::codon_ref_amino_acide)
        (*output) << '\t' << value.first << '/' << value.second;

    (*output) << "\tG+C\n";
}

void codon_count::add_to(ofstream *output, const map<string, int> &codon, const string &protName, int total, long gcNum, double gcpercent) {
    (*output) << protName;
    for (const auto &key: codon)
        (*output) << '\t' << key.second;
    (*output) << '\t' << gcNum;

    (*output) << "\nPercentage";
    for (const auto &key: codon)
        (*output) << '\t' << (((double)key.second / (double)total) * 100);
    (*output) << '\t' << gcpercent;

    (*output) << "\nRSCU";
    for (const auto &key: codon)
        (*output) << '\t' << protein::rscu(key.first, codon);
    (*output) << "\n\n";
}

int codon_count::check_args(const fs::path &input, const fs::path &output) {
    if (!exists(input)) {
        cout << "Cant't find directory at path : " << input << endl;
        return -1;
    }

    if (!directory::create_directories(output)) {
        cout << "Can't create or find directory at path : " << output << endl;
        return -1;
    }

    return 0;
}

int codon_count::main(const std::filesystem::path &input, const std::filesystem::path &output) {

    if (check_args(input, output) == -1) return -1;

    map<string, int> codon, total_codon;
    init_codon(codon);
    init_codon(total_codon);

    directory::to_fastaline(input);

    for (const auto &currentFile: fs::directory_iterator(input)) {
        if (!file::is_fastaline(currentFile)) continue;

        string fileName = currentFile.path().stem();
        ifstream *inputFile = file::read_open(currentFile.path());

        ofstream *outputFile = file::write_open((output / (fileName + ".tsv")), ios::trunc);
        init_output(outputFile);

        string lineRead, prot_name;
        int all_total(0);
        long gcTotal{0};
        size_t totalSize{0};
        while (getline((*inputFile), lineRead)) {
            if (lineRead.at(0) == '>') prot_name = lineRead;
            else {
                int total(0);
                for (unsigned long i = 0; i < lineRead.size(); i += 3) {
                    string sub = lineRead.substr(i, 3);
                    if (codon.find(sub) != codon.end()) {
                        codon[sub] += 1;
                        total_codon[sub] += 1;
                        total++;
                        all_total++;
                    } else {
                        cout << "Contig : " << prot_name << ", codon : " << sub << " unknown\n";
                    }
                }

                long gcCount = count_if(lineRead.begin(), lineRead.end(), [](char c) { return c == 'G' || c == 'C'; });
                double gcPourcent = ((double)gcCount / (double)lineRead.size()) * 100.0;

                add_to(outputFile, codon, prot_name, total, gcCount, gcPourcent);

                gcTotal += gcCount;
                totalSize += lineRead.size();

                reset_codon(codon);
            }
        }

        file::write_close(outputFile);
        file::read_close(inputFile);

        ofstream *total_output = file::write_open(output.string().append("/" + fileName + "-total.tsv"), ios::trunc);
        init_output(total_output);

        double totalPercentage = ((double)gcTotal / (double)totalSize) * 100.0;
        add_to(total_output, total_codon, fileName, all_total, gcTotal, totalPercentage);

        file::write_close(total_output);
        reset_codon(total_codon);
    }

    for (const auto &path: fs::directory_iterator(input)) {
        if (file::is_fastaline(path))
            remove(path);
    }

    return EXIT_SUCCESS;
}

void codon_count::init_codon(map<string, int> &codon) {
    for (auto &value: protein::codon_list)
        codon[value] = 0;
}

void codon_count::reset_codon(map<std::string, int> &codon) {
    for (auto &value: codon)
        value.second = 0;
}