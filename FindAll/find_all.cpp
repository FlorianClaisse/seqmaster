//
// Created by Florian Claisse on 02/05/2023.
//

#include <fstream>
#include "../Foundation/Headers/fasta.h"
#include "../Foundation/Headers/directory.h"
#include "find_all.h"

bool is_fasta_file(const std::filesystem::path &filePath) {
    return is_regular_file(filePath) && directory::have_extension(filePath, "fasta");
}

bool is_fastaline_file(const std::filesystem::path &filePath) {
    return is_regular_file(filePath) && directory::have_extension(filePath, "fastaline");
}

int find_all::start(const program_option::FindAll &options) {
    // convert inputA to fastaline
    if (!is_fasta_file(options.inputA)) {
        std::cout << "Le fichier : " << options.inputA <<  "n'est pas un fichier fasta." << std::endl;
        return EXIT_FAILURE;
    }
    if (fasta::to_fasta_line(options.inputA) != EXIT_SUCCESS) return EXIT_FAILURE;

    if (!std::filesystem::is_directory(options.inputB)) {
        std::cout << "Path : " << options.inputB << "n'est pas un dossier" << std::endl;
        return EXIT_FAILURE;
    }

    // Convert directory fasta to fastaline
    for (const auto &currentFile : std::filesystem::directory_iterator(options.inputB)) {
        if (!is_fasta_file(currentFile)) continue;
        if (fasta::to_fasta_line(currentFile) != EXIT_SUCCESS) return EXIT_FAILURE;
    }

    // Stocker les contigs du fichier de test dans un tableau.
    std::string inputPath(directory::removeExtension(options.inputA).append(".fastaline"));
    std::ifstream inputFile;
    inputFile.open(inputPath);

    std::string name;
    std::string value;
    std::vector<std::string> contigs_to_test_name;
    std::vector<std::string> contigs_to_test_value;
    while(!inputFile.eof()) {
        getline(inputFile, name);
        getline(inputFile, value);
        contigs_to_test_name.push_back(name);
        contigs_to_test_value.push_back(value);
    }
    inputFile.close();

    std::ofstream outputFile;
    outputFile.open(options.output.string().append("/output.txt"), std::ios::trunc);
    outputFile << "Filename\t";
    for (const auto &contig_name : contigs_to_test_name) {
        outputFile << contig_name << "\t";
    }
    outputFile << "\n";

    for (const auto &file : std::filesystem::directory_iterator(options.inputB)) {
        if (!is_fastaline_file(file)) continue;
        std::map<std::string, bool> results = fasta::find_contigs(file.path(), contigs_to_test_value);
        outputFile << file.path().filename() << "\t";
        for (const auto &contig_value : contigs_to_test_value) {
            outputFile << (results[contig_value] ? "V" : "X") << "\t";
        }
        outputFile << "\n";
    }

    outputFile.close();

    return EXIT_SUCCESS;
}