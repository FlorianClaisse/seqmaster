//
// Created by Florian Claisse on 23/05/2023.
//

#include <iostream>
#include <string>
#include <filesystem>

#include <sharg/all.hpp>

#include "include/parser.h"
#include "include/termcolor.hpp"

#include "../FindAll/include/find_all.hpp"
#include "../CodonCount/include/condo_count.h"
#include "../ContigDiff/include/contig_diff.h"
#include "../GeneMut/include/gene_mut.hpp"
#include "../ContigOrdered/include/contigOrdered.hpp"

#define APP_NAME "Contig"

#define FIND_ALL "findall"
#define CODON_COUNT "codoncount"
#define CONTIG_DIFF "contigdiff"
#define GENE_MUT "genemut"
#define CONTIG_ORDERED "contigordered"

using namespace std;
namespace fs = std::filesystem;

void initialise_parser(sharg::parser &parser) {
    parser.info.author = "Florian Claisse";
    parser.info.version = "0.0.1";
}

// =====================================================================================================================
// gene mut
// =====================================================================================================================

int run_contig_ordered(sharg::parser &parser) {
    initialise_parser(parser);
    fs::path input, output;

    parser.add_option(input, sharg::config{.short_id = 'i',
            .long_id = "input",
            .description = "Path to the input folder",
            .required = true});

    parser.add_option(output, sharg::config{.short_id = 'o',
            .long_id = "output",
            .description = "Path to the output folder",
            .required = true});

    try {
        parser.parse();
        return contig_ordered::main(input, output);
    } catch (sharg::parser_error const &ext) {
        std::cerr << termcolor::red << termcolor::bold
                  << "You managed to find" << CONTIG_ORDERED << ", it was already difficult, wasn't it?.\n"
                  << termcolor::cyan
                  << "Look at the help:\n"
                  << termcolor::reset
                  << ext.what()
                  << "\n";
        return -1;
    }
}

// =====================================================================================================================
// gene mut
// =====================================================================================================================

int run_gene_mut(sharg::parser &parser) {
    initialise_parser(parser);
    fs::path input, groupPath, output;

    parser.add_option(input, sharg::config{.short_id = 'i',
                                           .long_id = "input",
                                           .description = "Path to the input file (CSV file)",
                                           .required = true});

    parser.add_option(groupPath, sharg::config{.short_id = 'g',
            .long_id = "group",
            .description = "Path to the groups file (CSV file)",
            .required = true});

    parser.add_option(output, sharg::config{.short_id = 'o',
            .long_id = "output",
            .description = "Path to the output directory",
            .required = true});

    try {
        parser.parse();
        return gene_mut::main(input, groupPath, output);
    } catch (sharg::parser_error const &ext) {
        std::cerr << termcolor::red << termcolor::bold
                  << "You managed to find" << GENE_MUT << ", it was already difficult, wasn't it?.\n"
                  << termcolor::cyan
                  << "Look at the help:\n"
                  << termcolor::reset
                  << ext.what()
                  << "\n";
        return -1;
    }
}

// =====================================================================================================================
// contig diff
// =====================================================================================================================
struct contigdiff_argument {
    fs::path inputA;
    fs::path inputB;
    fs::path output;
    string type;
    int accept{100};
    int threads{4};
    unsigned long min_size{2};
};

void init_contig_diff(sharg::parser &parser, contigdiff_argument &args) {
    parser.add_option(args.inputA,
                      sharg::config{.short_id = 'A',
                                    .long_id = "inputA",
                                    .description = "Path to folder A.",
                                    .required = true});

    parser.add_option(args.inputB,
                      sharg::config{.short_id = 'B',
                                    .long_id = "inputB",
                                    .description = "Path to folder B.",
                                    .required = true});

    parser.add_option(args.output,
                      sharg::config{.short_id = 'o',
                                    .long_id = "output",
                                    .description = "Path to output folder.",
                                    .required = true});

    parser.add_option(args.type,
                      sharg::config{.short_id = 't',
                                    .long_id = "type",
                                    .description = "The file type (nucl/prot).",
                                    .required = true});

    parser.add_option(args.accept,
                      sharg::config{.long_id = "accept",
                                    .description = "Allows you to specify the minimum percentage to accept a contig as recognized."});

    parser.add_option(args.threads,
                      sharg::config{.long_id = "threads",
                                    .description = "Sets the number of threads the program can use"});

    parser.add_option(args.min_size,
                      sharg::config{.long_id = "minsize",
                                    .description = "Set the minimum size of sequences to find"});
}

int run_contig_diff(sharg::parser &parser) {
    initialise_parser(parser);
    contigdiff_argument args{};
    init_contig_diff(parser, args);

    try {
        parser.parse();
        return contig_diff::main(args.inputA, args.inputB, args.output, args.type, args.accept, args.threads, args.min_size);
    } catch (sharg::parser_error const &ext) {
        std::cerr << termcolor::red << termcolor::bold
                  << "You managed to find" << CONTIG_DIFF << ", it was already difficult, wasn't it?.\n"
                  << termcolor::cyan
                  << "Look at the help:\n"
                  << termcolor::reset
                  << ext.what()
                  << "\n";
        return -1;
    }
}


// =====================================================================================================================
// codon count
// =====================================================================================================================
struct codoncount_arguments {
    fs::path input;
    fs::path output;
};

void init_codon_count(sharg::parser &parser, codoncount_arguments &args) {
    parser.add_option(args.input,
                      sharg::config{.short_id = 'i',
                                    .long_id = "input",
                                    .description = "Path to the folder that contains the contig files where you have to count the codons.",
                                    .required = true});

    parser.add_option(args.output,
                      sharg::config{.short_id = 'o',
                                    .long_id = "output",
                                    .description = "Path to output folder.",
                                    .required = true});
}

int run_codon_count(sharg::parser &parser) {
    initialise_parser(parser);
    codoncount_arguments args{};
    init_codon_count(parser, args);

    try {
        parser.parse();
        return codon_count::main(args.input, args.output);
    } catch (sharg::parser_error const &ext) {
        std::cerr << termcolor::red << termcolor::bold
                  << "You managed to find" << CODON_COUNT << ", it was already difficult, wasn't it?.\n"
                  << termcolor::cyan
                  << "Look at the help:\n"
                  << termcolor::reset
                  << ext.what()
                  << "\n";
        return -1;
    }
}
// =====================================================================================================================
// find all
// =====================================================================================================================

struct findall_arguments {
    fs::path inputA;
    fs::path inputB;
    fs::path output;
    string type;
    int accept{100};
    int threads{4};
};

void init_find_all(sharg::parser &parser, findall_arguments &args) {
    parser.add_option(args.inputA,
                      sharg::config{.short_id = 'A',
                                    .long_id = "inputA",
                                    .description = "Path to the file that contains the contigs to find.",
                                    .required = true});

    parser.add_option(args.inputB,
                      sharg::config{.short_id = 'B',
                                    .long_id = "inputB",
                                    .description = "Path to the folder that contains the files where you have to find the contigs.",
                                    .required = true});

    parser.add_option(args.output,
                      sharg::config{.short_id = 'o',
                                    .long_id = "output",
                                    .description = "Path to the folder that will contain the output file(s).",
                                    .required = true});

    parser.add_option(args.type,
                      sharg::config{.short_id = 't',
                                    .long_id = "type",
                                    .description = "The file type (nucl/prot).",
                                    .required = true});

    parser.add_option(args.accept,
                      sharg::config{.long_id = "accept",
                                    .description = "Allows you to specify the minimum percentage to accept a contig as recognized."});

    parser.add_option(args.threads,
                      sharg::config{.long_id = "threads",
                                    .description = "Sets the number of threads the program can use."});

}

int run_find_all(sharg::parser &parser) {
    initialise_parser(parser);
    findall_arguments args{};
    init_find_all(parser, args);

    try {
        parser.parse();
        return find_all::main(args.inputA, args.inputB, args.output, args.type, args.accept, args.threads); // TODO: Change declaration
    } catch (sharg::parser_error const &ext) {
        std::cerr << termcolor::red << termcolor::bold
                  << "You managed to find" << FIND_ALL << ", it was already difficult, wasn't it?.\n"
                  << termcolor::cyan
                  << "Look at the help:\n"
                  << termcolor::reset
                  << ext.what()
                  << "\n";
        return -1;
    }
}

// =====================================================================================================================
// main
// =====================================================================================================================
int parser::parse(int argc, const char **argv) {
    sharg::parser top_level_parser{APP_NAME, argc, argv, sharg::update_notifications::on, {FIND_ALL, CODON_COUNT, CONTIG_DIFF, GENE_MUT, CONTIG_ORDERED}};
    initialise_parser(top_level_parser);

    try {
        top_level_parser.parse();
    } catch (sharg::parser_error const &ext) {
        std::cerr << termcolor::red << termcolor::bold
                  << "OHHH no you don't know how to write the command line for this program.\n"
                  << termcolor::cyan
                  << "Look at the help:\n"
                  << termcolor::reset
                  << ext.what()
                  << "\n";

        return -1;
    }

    sharg::parser &sub_parser = top_level_parser.get_sub_parser();

    if (sub_parser.info.app_name == string_view{string{APP_NAME} + "-" + FIND_ALL}) {
        return run_find_all(sub_parser);
    } else if (sub_parser.info.app_name == string_view{string{APP_NAME} + "-" + CODON_COUNT}) {
        return run_codon_count(sub_parser);
    } else if (sub_parser.info.app_name == string_view{string{APP_NAME} + "-" + CONTIG_DIFF}) {
        return run_contig_diff(sub_parser);
    } else if (sub_parser.info.app_name == string_view{string{APP_NAME} + "-" + GENE_MUT}) {
        return run_gene_mut(sub_parser);
    } else if (sub_parser.info.app_name == string_view{string{APP_NAME} + "-" + CONTIG_ORDERED}) {
        return run_contig_ordered(sub_parser);
    } else return -1;
}