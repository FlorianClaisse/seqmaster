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
    int check_args(const fs::path &input, const fs::path &groupPath, const fs::path &output);
    void decode_group_file(const fs::path &groupPath, unordered_map<string, string> &genes);
}

int gene_mut::check_args(const fs::path &input, const fs::path &groupPath, const fs::path &output) {
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

int gene_mut::main(const fs::path &input, const fs::path &groupPath, const fs::path &output) {
    if (check_args(input, groupPath, output) == -1) return -1;
    // [gene_name, group_name]
    unordered_map<string, string> genes;
    decode_group_file(groupPath, genes);

    // [group_name, [gene_name]]
    unordered_map<string, unordered_set<string>> groups;
    for (const auto &value: genes)
        groups[value.second].insert(value.first);

    CSVReader reader{input.string()};

    int pos{reader.begin()->operator[](POS_REF).get<int>()};
    string gene{reader.begin()->operator[](ENTRY).get<string>()};
    string group{genes.find(gene)->second};
    size_t groupSize{groups.find(group)->second.size()};
    for (const auto &row: reader) {
        if (row[POS_REF] == pos) { // same on differente gene

        }
    }

    return 0;
}

void gene_mut::decode_group_file(const fs::path &groupPath, unordered_map<string, string> &genes) {
    CSVReader reader{groupPath.string()};

    vector<string> colName{reader.get_col_names()};
    for (auto &row: reader) {
        for (int i = 0; i < row.size(); i++)
            genes[row[i].get()] = colName[i];
    }
}