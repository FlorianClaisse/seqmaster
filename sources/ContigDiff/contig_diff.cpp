//
// Created by Florian Claisse on 10/05/2023.
//

#include "include/contig_diff.h"
#include "include/search.hpp"

#include "../Utils/include/termcolor.hpp"
#include "../Utils/include/directory.h"
#include "../Utils/include/file.h"

#define PROT "prot"
#define NUCL "nucl"

#define FIND_COMMON "find_common.fasta"

using namespace std;
namespace fs = std::filesystem;

namespace contig_diff {
    int error_message(const string &message);
    int check_options(const fs::path &inputA, const fs::path &inputB, const fs::path &output, const string &type, int accept, int threads);
}

int contig_diff::error_message(const string &message) {
    cout << termcolor::red << termcolor::bold
         << message
         << termcolor::reset;
    return -1;
}

int contig_diff::check_options(const fs::path &inputA, const fs::path &inputB, const fs::path &output, const string &type, int accept, int threads) {
    if (threads <= 0) return error_message("You want to run the program on a number of threads less than or equal to 0, go try again ;).\n");
    if (type != PROT && type != NUCL) return error_message("For the type, the value must be \"prot\" or \"nucl\" did you read the docs ?.\n");
    if (accept < 0 || accept > 100) return error_message("You just gave an acceptance percentage lower than 0 or higher than 100. That's not very clever.\n");
    if (!is_directory(inputA)) return error_message("Input A is not a folder or does not exist. Try again the next one is the right one.\n");
    if (!is_directory(inputB)) return error_message("Input B is not a folder or does not exist. Try again the next one is the right one.\n");
    if(!directory::create_directories(output)) return error_message("Unable to find/create the output folder are you sure you have given a valid path.\n");

    return 0;
}

template<typename traits_t, typename config_t>
int start_search(const contig_diff::param &options, const config_t &config) {

    seqan3::sequence_file_output f_out{options.output / FIND_COMMON};

    cout << "Start all common in A.\n";
    contig_diff::Search<traits_t, config_t> search(options, config);

    search.search_common(options.inputA);

    return 0;
}

int contig_diff::main(const fs::path &inputA, const fs::path &inputB, const fs::path &output, const string &type, int accept, int threads) {

    if (check_options(inputA, inputB, output, type, accept, threads) != 0) return -1;

    int nb_fasta = directory::count_fasta_file(inputA);
    if (nb_fasta < 2) {
        cout << "Il faut au minimum deux fichiers de type fasta dans le dossier A.\n";
        return -1;
    }
    nb_fasta = directory::count_fasta_file(inputB);
    if (nb_fasta < 1) {
        cout << "Il faut au minimum 1 fichier de type fasta dans le dossier B.\n";
        return -1;
    }

    param options = {inputA, inputB, output, (type == NUCL), (100 - accept), threads};

    if (options.nucl) {
        seqan3::configuration const config = seqan3::align_cfg::method_global{
                seqan3::align_cfg::free_end_gaps_sequence1_leading{false},
                seqan3::align_cfg::free_end_gaps_sequence2_leading{true},
                seqan3::align_cfg::free_end_gaps_sequence1_trailing{false},
                seqan3::align_cfg::free_end_gaps_sequence2_trailing{true}}
                                             | seqan3::align_cfg::scoring_scheme{seqan3::nucleotide_scoring_scheme{
                seqan3::match_score{0},
                seqan3::mismatch_score{-1}}}
        | seqan3::align_cfg::parallel{static_cast<uint32_t>(options.threads)};

        using traits_t = seqan3::sequence_file_input_default_traits_dna;

        start_search<traits_t>(options, config);
    } else {
        seqan3::configuration const config = seqan3::align_cfg::method_global{
                seqan3::align_cfg::free_end_gaps_sequence1_leading{false},
                seqan3::align_cfg::free_end_gaps_sequence2_leading{true},
                seqan3::align_cfg::free_end_gaps_sequence1_trailing{false},
                seqan3::align_cfg::free_end_gaps_sequence2_trailing{true}}
                                             | seqan3::align_cfg::scoring_scheme{seqan3::aminoacid_scoring_scheme{
                seqan3::match_score{0},
                seqan3::mismatch_score{-1}}}
        | seqan3::align_cfg::parallel{static_cast<uint32_t>(options.threads)};

        using traits_t = seqan3::sequence_file_input_default_traits_aa;

        start_search<traits_t>(options, config);
    }

    return 0;

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

/*void contig_diff::generate_output(const fs::path &directory_path, const map<std::string, contig_diff::Common> &common) {
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
}*/


