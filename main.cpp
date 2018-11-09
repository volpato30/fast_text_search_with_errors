#include <iostream>
#include <fstream>
#include <algorithm>
#include <chrono>

#include <assert.h>
#include "hamming1_search.h"
#include "string_split.h"


void test_hamming1_search() {
    auto pattern = std::vector<std::string> {"T", "C", "G", "T"};
    auto text = std::vector<std::string> {"T",
                                          "T",
                                          "T",
                                          "A",
                                          "C",
                                          "G",
                                          "T",
                                          "A",
                                          "A",
                                          "A",
                                          "C",
                                          "T",
                                          "A",
                                          "A",
                                          "A",
                                          "C",
                                          "T",
                                          "G",
                                          "T",
                                          "A",
                                          "A"};
    auto match_result = hamming1_search(pattern, text);
    assert(!match_result.exactMatch);
    assert(match_result.startIndex.size() == 1);
    assert(match_result.startIndex[0] == 3);

    pattern = std::vector<std::string> {"A", "A", "C", "T"};
    match_result = hamming1_search(pattern, text);
    assert(match_result.exactMatch);
}

void readIdentifiedFile(const char * filename, std::vector<std::vector<std::string>> &identifiedPeptides) {
    identifiedPeptides.clear();
    std::ifstream file(filename);
    std::string s;
    while (getline(file, s)) {
        s.erase(std::remove(s.begin(), s.end(), '\n'), s.end());
        identifiedPeptides.push_back(split(s, ','));
    }
    file.close();
}

std::string match(const std::vector<std::string> &pattern, const std::vector<std::vector<std::string>> &identifiedPeptides) {
    matchLocations ml;
    std::string output_string;
    unsigned long m = pattern.size();
    std::vector<std::string> wildTypePeptideVector = std::vector<std::string>();

    for (const std::vector<std::string> &wildTypePeptide : identifiedPeptides) {
        ml = hamming1_search(pattern, wildTypePeptide);
        if (ml.exactMatch) {
            continue;
        } else if (ml.startIndex.empty()) {
            continue;
        } else {
            for (const unsigned long &start_location : ml.startIndex) {
                wildTypePeptideVector.push_back(join_subvector(wildTypePeptide, ',', start_location, m));
            }
        }
    }
    if (!wildTypePeptideVector.empty()) {
        join(wildTypePeptideVector, ';', output_string);
        output_string = join(pattern, ',') + '\t' + output_string + '\n';
    }

    return output_string;
}

int main() {
    // testing
    test_hamming1_search();
    std::cout << "passed test!" << std::endl;
    char identified_filename[] = "/home/rui/CLionProjects/wild-type-match/identified_features.cc.txt";
    char output_filename[] = "/home/rui/CLionProjects/wild-type-match/wildtype_matched.peptides.txt";
    char denovo_filename[] = "/home/rui/CLionProjects/wild-type-match/denovo_peptides.cc.txt";

    // read in identified peptides;
    auto identifiedPeptides = std::vector<std::vector<std::string>> ();

    readIdentifiedFile(identified_filename, identifiedPeptides);
    std::cout << identifiedPeptides.size() << " identified peptides" << std::endl;
    std::cout << identifiedPeptides[1].size() << std::endl;

    // match wild-type and write to output
    int counter = 0;

    std::ifstream denovo_file(denovo_filename);
    std::string denovo_peptides, line_string;
    std::vector<std::string> pattern;

    std::ofstream outputFile(output_filename);
    auto t1 = std::chrono::high_resolution_clock::now();
    auto t2 = t1;
    while (getline(denovo_file, denovo_peptides)) {
        // strip \n
        counter++;
        if (counter % 1000 == 0) {
            std::cout << "processed " << counter << "peptides" << std::endl;
            t2 = std::chrono::high_resolution_clock::now();
            std::cout << "takes " << std::chrono::duration_cast<std::chrono::seconds>(t2-t1).count() << " seconds" << std::endl;
            t1 = t2;
        }

        denovo_peptides.erase(std::remove(denovo_peptides.begin(), denovo_peptides.end(), '\n'), denovo_peptides.end());
        pattern.clear();
        pattern = split(denovo_peptides, ',');
        line_string = match(pattern, identifiedPeptides);
        if (!line_string.empty()) {
            std::cout << "write: " << line_string << std::endl;
            if (outputFile.is_open()) {
                outputFile << line_string;
                outputFile.flush();
            } else {
                std::cerr << "not writing!!" << std::endl;
            }
        }
    }
    outputFile.close();
    denovo_file.close();
    return 0;
}