//
// Created by Florian Claisse on 02/05/2023.
//

#include <fstream>
#include "../Foundation/Headers/fasta.h"
#include "../Foundation/Headers/directory.h"
#include "find_all.h"

using namespace std;
namespace fs = std::filesystem;

int find_all::start(const program_option::FindAll &options) {
    // convert inputA to fastaline
    if (!fasta::is_fasta_file(options.inputA)) {
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
        if (!fasta::is_fasta_file(currentFile)) continue;
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
    remove(fs::path(inputPath));

    ofstream outputFile, currentOutputResult;
    outputFile.open(options.output.string().append("/output.txt"), ios::trunc);
    outputFile << "Filename\t\n";

    for (const auto &file : filesystem::directory_iterator(options.inputB)) {
        if (!fasta::is_fastaline_file(file)) continue;

        map<string, bool> results;
        if (options.accept == 100) results = fasta::find_contigs(file.path(), contigs_to_test_value);
        else results = fasta::find_contigs(file.path(), contigs_to_test_value, options.accept);

        string fileNameWithoutExtension(directory::fileNameWithoutExtension(file.path()));
        currentOutputResult.open(options.output.string().append("/" + fileNameWithoutExtension + "-result.fasta"), ios::trunc),

        outputFile << directory::fileNameWithoutExtension(file.path()) << "\t";
        for (const auto &contig_value : contigs_to_test_value) {
            if (results[contig_value]) {
                int index = &contig_value - &contigs_to_test_value[0];
                outputFile << (contigs_to_test_name[index]) << "\t";
                currentOutputResult << (contigs_to_test_name[index]) << endl << contig_value << endl;
            }
        }
        outputFile << "\n";
        currentOutputResult.close();
        remove(file);
    }

    outputFile.close();

    return EXIT_SUCCESS;
}