//
// Created by Florian Claisse on 02/05/2023.
//

#include <fstream>
#include <filesystem>
#include <omp.h>
#include "../Foundation/include/fasta.h"
#include "../Foundation/include/directory.h"
#include "find_all.h"

using namespace std;
namespace fs = std::filesystem;

int check_options(const program_option::FindAll &options) {
    if (!fasta::is_fasta_file(options.inputA)) {
        cout << "Le fichier : " << options.inputA <<  "n'est pas un fichier fasta." << endl;
        return EXIT_FAILURE;
    }

    if (!fs::is_directory(options.inputB)) {
        cout << "Path : " << options.inputB << "n'est pas un dossier" << endl;
        return EXIT_FAILURE;
    }

    if (!fs::is_directory(options.output)) {
        cout << "Path : " << options.output << "n'est pas un dossier" << endl;
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}

int find_all::start(const program_option::FindAll &options) {
    if (check_options(options) != EXIT_SUCCESS) return EXIT_FAILURE;

    if (fasta::to_fasta_line(options.inputA) != EXIT_SUCCESS) return EXIT_FAILURE;

    // Convert directory fasta to fastaline
    for (const auto &currentFile : fs::directory_iterator(options.inputB)) {
        if (!fasta::is_fasta_file(currentFile)) continue;
        if (fasta::to_fasta_line(currentFile) != EXIT_SUCCESS) return EXIT_FAILURE;
    }

    // Stocker les contigs du fichier de test dans un tableau.
    string inputPath(directory::removeExtension(options.inputA).append(".fastaline"));
    ifstream inputFile;
    inputFile.open(inputPath);

    string name;
    string value;
    map<string, string> contigs;
    while(!inputFile.eof()) {
        getline(inputFile, name);
        getline(inputFile, value);
        contigs[name] = value;
    }
    for (const auto &contig_value: contigs) {
        cout << "Value : " << contig_value.first << ", name : " << contig_value.second << endl << endl;
    }
    inputFile.close();

    ofstream outputFile;
    outputFile.open(options.output.string().append("/output.txt"), ios::trunc);
    outputFile << "Filename\t\n";

    for (const auto &file : filesystem::directory_iterator(options.inputB)) {

        if (!fasta::is_fastaline_file(file)) continue;

        string fileNameWithoutExtension(directory::fileNameWithoutExtension(file.path()));
        ofstream currentOutputResult;
        currentOutputResult.open(options.output.string().append("/" + fileNameWithoutExtension + "-result.fasta"), ios::trunc);

        outputFile << directory::fileNameWithoutExtension(file.path()) << "\t";

        if (options.accept == 100) {
            fasta::find_contig(file.path(), contigs, options.nucl, [&outputFile, &currentOutputResult](const string &nameA, const string &nameB, const string &value) -> void {
                outputFile << nameA << "\t";
                currentOutputResult << nameB << " -> " << nameA << endl << value << endl;
            });
        }
        else {
            fasta::find_contigs(file.path(), contigs, (100 - options.accept), options.nucl, [&outputFile, &currentOutputResult](const string &nameA, const string &nameB, const string &value, double percentage) -> void {
                outputFile << nameA << "\t";
                currentOutputResult << nameB << " -> " << nameA << " -> " << (100.0 - percentage) << "%" << endl << value << endl;
            });
        }

        outputFile << "\n";
        currentOutputResult.close();
        remove(file);
    }

    outputFile.close();

    return EXIT_SUCCESS;
}