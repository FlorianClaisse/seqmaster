//
// Created by Florian Claisse on 02/05/2023.
//

#include <fstream>
#include <filesystem>
#include <iostream>
#include <omp.h>
#include "../Utils/include/fasta.h"
#include "../Utils/include/directory.h"
#include "include/find_all.h"

using namespace std;
namespace fs = std::filesystem;

// Format : >ContigNameB -> ContigNameA -> Percentage
int generate_output(const find_all::option &options) {
    ofstream outputFile(options.output.string().append("/output.txt"), ios::trunc);
    outputFile << "Filename\n";
    string lineRead;
    for (const auto &currentPath : fs::directory_iterator(options.output)) {
        if (!directory::is_result_file(currentPath)) continue;
        outputFile << currentPath.path().filename() << "\t";
        ifstream currentFile(currentPath.path());
        while(getline(currentFile, lineRead)) {
            if (lineRead.empty()) continue;
            if (lineRead.at(0) == '>') {
                size_t pos, second_pos;
                if ((pos = lineRead.find("->")) != string::npos) {
                    outputFile << lineRead.substr(pos + 3) << "\t";
                } else return EXIT_FAILURE;
            }
        }
        outputFile << "\n";
    }

    outputFile.close();

    return EXIT_SUCCESS;
}

int find_all::main(const find_all::option &options) {

    // Convert all input to fastaline file.
    if (fasta::to_fastaline(options.inputA) == EXIT_FAILURE) return EXIT_FAILURE;
    if (fasta::directory_to_fastaline(options.inputB) == EXIT_FAILURE) return EXIT_FAILURE;

    ifstream *inputFile = directory::read_open(fs::path(options.inputA).replace_extension("fastaline"));
    // Stocker les contigs du fichier de test dans un tableau.
    string name, value;
    map<string, string> contigs;
    while(!inputFile->eof()) {
        getline((*inputFile), name);
        getline((*inputFile), value);
        contigs[name] = value;
    }

    directory::read_close(inputFile);

    long nb_files = distance(fs::directory_iterator(options.inputB), fs::directory_iterator{}); // C++17
    cout << "Nb files : " << nb_files << endl;

#pragma omp parallel for default(none) shared(nb_files, contigs, options) num_threads(options.threads)
    for (long i = 0; i < nb_files; i++) {
        fs::directory_entry file = *next(fs::directory_iterator(options.inputB), i);
        if (!directory::is_fastaline_file(file)) continue;
        if (directory::is_result_file(file)) continue;

        printf("File : %s, Threads = %d, i = %ld\n", file.path().c_str(), omp_get_thread_num(), i);

        string fileNameWithoutExtension(file.path().stem());
        ofstream *currentOutputResult = directory::write_open(options.output.string().append("/" + fileNameWithoutExtension + "-result.fasta"), ios::trunc);

        if (options.accept == 100) {
            fasta::find_contig(file.path(), contigs, options.nucl, [&currentOutputResult](const string &nameA, const string &nameB, const string &value) -> void {
                (*currentOutputResult) << nameB << " -> " << nameA << endl << value << endl;
            });
        }
        else {
            fasta::find_contigs(file.path(), contigs, (100 - options.accept), options.nucl, [&currentOutputResult](const string &nameA, const string &nameB, const string &value, double percentage) -> void {
                (*currentOutputResult) << nameB << " -> " << nameA << " -> " << (100.0 - percentage) << "%" << endl << value << endl;
            });
        }

        directory::write_close(currentOutputResult);
    }

    return generate_output(options);
}