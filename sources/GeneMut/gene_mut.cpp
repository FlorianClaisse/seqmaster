//
// Created by Florian Claisse on 14/06/2023.
//

#include <iostream>
#include <unordered_map>
#include <unordered_set>

#include "csv.hpp"
#include "../Utils/include/directory.h"

#include "include/gene_mut.hpp"

#define ENTRY "Entry"
#define LABEL "Label"
#define POS_REF "Position reference"

using namespace std;
using namespace csv;
namespace fs = std::filesystem;

namespace gene_mut {
    int check_args(const fs::path &input, const fs::path &groupPath, const fs::path &output) {
        if (!exists(groupPath)) {
            cout << "Cannot find file at path : " << groupPath << endl;
            return -1;
        }

        if (!exists(input)) {
            cout << "Cannot find file at path : " << input << endl;
            return -1;
        }

        if (!directory::create_directories(output)) {
            cout << "Failed to find or create output directory at path : " << output << endl;
            return -1;
        }

        return 0;
    }

    void decode_group_file(const fs::path &groupPath, unordered_map<string, string> &genes) {
        CSVReader reader{groupPath.string()};

        vector<string> colName{reader.get_col_names()};
        for (auto &row: reader) {
            for (int i = 0; i < row.size(); i++)
                genes[row[i].get()] = colName[i];
        }
    }

    string groupName(const string &gene, const unordered_map<string, string> &genes) {
        auto fValue = genes.find(gene);

        if (fValue == genes.end()) throw invalid_argument("Can't retrieve : " + gene + " in gene groups");
        return fValue->second;
    }

    void reset_rows(const CSVRow &row, unordered_map<string, vector<CSVRow>> &rows, int &pos, const unordered_map<string, string> &genes) {
        string group_name = groupName(row[ENTRY].get<string>(), genes);

        rows.clear();
        rows[group_name].push_back(row);
        pos = row[POS_REF].get<int>();
    }

    void add_to_rows(const CSVRow &row, unordered_map<string, vector<CSVRow>> &rows, const unordered_map<string, string> &genes) {
        string group_name = groupName(row[ENTRY].get<string>(), genes);

        rows[group_name].push_back(row);
    }

    vector<string> row_to_array(const CSVRow &row) {
        vector<string> result;
        for (auto &field: row)
            result.push_back(field.get<string>());

            return result;
    }
}

int gene_mut::main(const fs::path &input, const fs::path &groupPath, const fs::path &output) {
    if (check_args(input, groupPath, output) == -1) return -1;
    // [gene_name, group_name]
    unordered_map<string, string> genes;
    decode_group_file(groupPath, genes);

    // [group_name, [gene_name]]
    unordered_map<string, vector<string>> groups;
    for (const auto &value: genes)
        groups[value.second].push_back(value.first);

    CSVReader reader{input.string()};

    // [group_name, output_file]
    unordered_map<string, CSVWriter<ofstream>*> outputs;
    for (const auto &group: groups) {
        auto *ou_t = new ofstream{output / (group.first + ".csv"), ios::trunc};
        outputs[group.first] = new CSVWriter<ofstream>{(*ou_t)};
        (*outputs[group.first]) << reader.get_col_names();
    }

    // [group_name, [rows]]
    unordered_map<string, vector<CSVRow>> currentRows;
    int currentPos{0};
    reset_rows((*reader.begin()), currentRows, currentPos, genes);
    for (const auto &row: reader) {
        if (row[POS_REF].get<int>() == currentPos)
            add_to_rows(row, currentRows, genes);
        else {
            for (const auto &value: currentRows) {
                auto fValue = groups.find(value.first);
                if (fValue == groups.end()) throw invalid_argument("Can't find " + value.first + " inside groups");

                if (value.second.size() == fValue->second.size()) {
                    auto sfValue = outputs.find(value.first);
                    if (sfValue == outputs.end()) throw invalid_argument("Can't find " + value.first + " inside outputs");

                    CSVWriter<ofstream> *ou_t{sfValue->second};
                    for (const auto &write_row: value.second) {
                        (*ou_t) << row_to_array(write_row);
                    }
                }
            }

            reset_rows(row, currentRows, currentPos, genes);
        }
    }

    for (const auto &value: outputs)
        delete value.second;

    return 0;
}