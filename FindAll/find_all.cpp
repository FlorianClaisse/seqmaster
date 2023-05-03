//
// Created by Florian Claisse on 02/05/2023.
//

#include <fstream>
#include "../Foundation/Headers/fasta.h"
#include "../Foundation/Headers/directory.h"
#include "find_all.h"

using namespace std;
namespace fs = std::filesystem;

bool is_fasta_file(const fs::path &filePath) {
    return is_regular_file(filePath) && directory::have_extension(filePath, "fasta");
}

bool is_fastaline_file(const fs::path &filePath) {
    return is_regular_file(filePath) && directory::have_extension(filePath, "fastaline");
}

int find_all::start(const program_option::FindAll &options) {
    // convert inputA to fastaline
    if (!is_fasta_file(options.inputA)) {
        cout << "Le fichier : " << options.inputA <<  "n'est pas un fichier fasta." << endl;
        return EXIT_FAILURE;
    }
    if (fasta::to_fasta_line(options.inputA) != EXIT_SUCCESS) return EXIT_FAILURE;

    if (!fs::is_directory(options.inputB)) {
        cout << "Path : " << options.inputB << "n'est pas un dossier" << endl;
        return EXIT_FAILURE;
    }

    // Convert directory fasta to fastaline
    for (const auto &currentFile : fs::directory_iterator(options.inputB)) {
        if (!is_fasta_file(currentFile)) continue;
        if (fasta::to_fasta_line(currentFile) != EXIT_SUCCESS) return EXIT_FAILURE;
    }

    // Stocker les contigs du fichier de test dans un tableau.
    string inputPath(directory::removeExtension(options.inputA).append(".fastaline"));
    ifstream inputFile;
    inputFile.open(inputPath);

    string name;
    string value;
    vector<string> contigs_to_test_name;
    vector<string> contigs_to_test_value;
    while(!inputFile.eof()) {
        getline(inputFile, name);
        getline(inputFile, value);
        contigs_to_test_name.push_back(name);
        contigs_to_test_value.push_back(value);
    }
    inputFile.close();

    ofstream outputFile;
    outputFile.open(options.output.string().append("/output.txt"), ios::trunc);
    outputFile << "Filename\t";
    for (const auto &contig_name : contigs_to_test_name) {
        outputFile << contig_name << "\t";
    }
    outputFile << "\n";

    for (const auto &file : filesystem::directory_iterator(options.inputB)) {
        if (!is_fastaline_file(file)) continue;
        map<string, bool> results = fasta::find_contigs(file.path(), contigs_to_test_value);
        outputFile << file.path().filename() << "\t";
        for (const auto &contig_value : contigs_to_test_value) {
            outputFile << (results[contig_value] ? "V" : "X") << "\t";
        }
        outputFile << "\n";
    }

    outputFile.close();

    return EXIT_SUCCESS;
}