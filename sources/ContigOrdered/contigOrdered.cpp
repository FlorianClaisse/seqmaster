//
// Created by Florian Claisse on 21/06/2023.
//

#include <iostream>
#include <string>
#include <vector>
#include <cstring>

#include "../Utils/include/file.h"
#include "../Utils/include/directory.h"

#include "include/contigOrdered.hpp"

using namespace std;
namespace fs = std::filesystem;

namespace contig_ordered {
    int check_argument(const fs::path &input, const fs::path &output) {
        if (!exists(input)) {
            cout << "Can't find directory at path : " << input << '\n';
            return -1;
        }
        if (!directory::create_directories(output)) {
            cout << "Can't create or find directory at path : " << output << '\n';
            return -1;
        }

        for (const auto &filePath: fs::directory_iterator(output))
            fs::remove_all(filePath);

        return 0;
    }

    string contig_number(int i) {
        if (i >= 1 && i < 10) return (">Contig_000" + to_string(i));
        if (i >= 10 && i < 100) return (">Contig_00" + to_string(i));
        if (i >= 100 && i < 1000) return (">Contig_0" + to_string(i));
        else return (">Contig_" + to_string(i));
    }
}

int contig_ordered::main(const std::filesystem::path &input, const std::filesystem::path &output) {
    if (check_argument(input, output) == -1) return -1;
    directory::clean_fastas(input);

    for (const auto &path: fs::directory_iterator{input}) {
        if (!file::is_fasta(path)) continue;

        ifstream in_f{path.path()};

        string line;
        vector<string> contigs{};
        while(getline(in_f, line)) {
            if (line.empty()) continue;
            if (line.at(0) == '>') continue;

            transform(line.begin(), line.end(), line.begin(), ::toupper);
            stringstream ss{line};

            string sub;
            while(getline(ss, sub, 'N')) {
                if (!sub.empty())
                    contigs.push_back(sub);
            }
        }
        in_f.close();

        ofstream ou_f{output / path.path().filename()};
        for (int i = 0; i < contigs.size(); i++)
            ou_f << contig_number(i) << '\n' << contigs[i] << '\n';
        ou_f.close();
    }

    return 0;
}