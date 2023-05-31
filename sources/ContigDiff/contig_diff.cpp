//
// Created by Florian Claisse on 10/05/2023.
//

#include <iostream>
#include <filesystem>
#include <fstream>
#include <map>
#include <vector>
#include <array>
#include <algorithm>

#include "../Utils/include/fasta.h"
#include "../Utils/include/directory.h"
#include "../Utils/include/file.h"

#include "include/contig_diff.h"
#include "include/cent_percent.h"
#include "include/other_percent.h"

// TODO: Optimisation: MultiThreads
// TODO: Amelioreation: Prendre en compte un pourcentage different de 100%

using namespace std;
namespace fs = std::filesystem;

namespace contig_diff {
    void generate_output(const fs::path &directory_path, const map<string, Common> &common);
}

int contig_diff::main(const std::filesystem::path &inputA, const std::filesystem::path &inputB, const std::filesystem::path &output,
                      const std::string &type, int accept, int threads) {
    /*cout << "Convert input A file's to fastaline.\n";
    fasta::directory_to_fastaline(options.inputA);

    cout << "Convert input B file's to fastaline.\n";
    fasta::directory_to_fastaline(options.inputB);

    int nb_fasta_file = directory::count_fasta_file(options.inputA);
    if (nb_fasta_file < 2) {
        cout << "Il faut au minimum deux fichiers de type fasta dans le dossier A.\n";
        return EXIT_FAILURE;
    }

    map<string, Common> common;

    if (options.accept == 100) common = contig_diff::cent_percent::main(options);
    //else common = contig_diff::other_percent::main(options);

    cout << "Start Generate output\n";
    generate_output(options.output, common);
    cout << "Finish\n";
    init_max_nb_error(5);

    // [sub_contig_value, Common]

    // Init common with 2 file
    search_with_2_file(two_first.first, two_first.second, common2, 5);

    // Compare first file with all files
    for (const auto &path: fs::directory_iterator(options.inputA)) {
        if (path == two_first.first || !directory::is_fastaline_file(path)) continue;
        search_with_2_file(path, two_first.first, common, 100);
        search_with_2_file(path, two_first.first, common2, 5);
    }

    cout << "Size1 : " << common.size() << ", Size2 : " << common2.size() << endl;
    cout << "Equal : " << (common == common2) << endl;

    generate_output(options.output, common);
    generate_output("/Users/florianclaisse/Documents/Test/contigdiff/", common2);

    // Check all possible common inside all file in A.
    check_on_directory(options.inputA, common);
    cout << "Finish\nAll common inside A size : " << common.size() << endl;
    cout << "Start check all common in A with B.\n";
    // Check all common in A with B;
    verif_on_directory(options.inputB, common);
    cout << "Finish\nAll common after check in B size : " << common.size() << endl;
    cout << "Start Generate output\n";

    cout << "Finish\n";

    return EXIT_SUCCESS;*/
}

/*pair<unsigned long, double> contig_diff::find_inside_file(ifstream &test_file, const string &word, int accept) { // Verifier
    unsigned long max_size(0), nb_errors;
    double percentage(0), current_percentage;

    string line_read;
    while(getline(test_file, line_read)) {
        if (line_read.at(0) == '>') continue;
        for (unsigned long offset = 0; offset < line_read.size(); offset++) {
            nb_errors = 0;
            for (unsigned long i = 0; i < word.size(); i++) {
                if (offset + i >= line_read.size()) {
                    for (unsigned long j = i; j < word.size(); j++) {
                        nb_errors++;
                        if (nb_errors <= max_nb_error[j] && j >= max_size) {
                            current_percentage = (100 * (double)nb_errors) / max_nb_error[j];
                            if (j == max_size && current_percentage < percentage) percentage = current_percentage;
                            else {
                                max_size = j;
                                percentage = current_percentage;
                            }
                        }
                    }
                    break;
                }
                if (line_read[i + offset] != word[i]) nb_errors++;
                if (nb_errors <= max_nb_error[i] && i >= max_size) {
                    current_percentage = (100 * (double)nb_errors) / max_nb_error[i];
                    if (i == max_size && current_percentage < percentage) percentage = current_percentage;
                    else {
                        max_size = i;
                        percentage = current_percentage;
                    }
                }
            }
        }
    }

    return {(max_size + 1), percentage};
}*/

void contig_diff::generate_output(const fs::path &directory_path, const map<std::string, contig_diff::Common> &common) {
    string outputPath = directory_path.string().append("/output.txt");
    ofstream *output = file::write_open(outputPath, ios::trunc);
    (*output) << "Filename\tContig name\tContig value\tPercentage\tFind in B\tA - B\n";

    for (const auto &value: common) {
        (*output) << value.second.filename << "\t"
                  << value.second.contigName << "\t"
                  << value.first << "\t"
                  << value.second.percentage << "\t"
                  << value.second.inputB << "\t"
                  << value.first.substr(value.second.inputB.size()) << "\t\n";
    }

    file::write_close(output);
}


