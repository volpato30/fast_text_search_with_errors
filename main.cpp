#include <iostream>
#include <fstream>
#include <algorithm>
#include <chrono>
#include <future>
#include <mutex>
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

std::string match(const std::vector<std::string> &pattern, const std::vector<std::vector<std::string>> &identifiedPeptides,
                  const std::function<std::string(std::string)> &map_func=string_identity_map) {
    matchLocations ml;
    std::string output_string;
    unsigned long m = pattern.size();
    std::vector<std::string> wildTypePeptideVector = std::vector<std::string>();

    for (const std::vector<std::string> &wildTypePeptide : identifiedPeptides) {
        ml = hamming1_search(pattern, wildTypePeptide, map_func);
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

void _process_one_line(std::string denovo_peptide, const std::vector<std::vector<std::string>> *identifiedPeptides,
        std::ofstream *outputFile, std::mutex &m) {
    denovo_peptide.erase(std::remove(denovo_peptide.begin(), denovo_peptide.end(), '\n'), denovo_peptide.end());
    std::vector<std::string> pattern = split(denovo_peptide, ',');
    std::string line_string = match(pattern, *identifiedPeptides, I_to_L_map);

    if (!line_string.empty()) {
        std::cout << "write: " << line_string << std::endl;
        m.lock();
        if (outputFile->is_open()) {
            *outputFile << line_string;
            outputFile->flush();
        } else {
            std::cerr << "not writing!!" << std::endl;
        }
        m.unlock();
    }
}

void parallel_match(std::vector<std::string>::iterator beg, std::vector<std::string>::iterator end,
        const std::vector<std::vector<std::string>> *identifiedPeptides, std::ofstream *outputFile, std::mutex &m) {
    std::cout << "start a new thread\n";
    auto len = end - beg;
    if (len < 2000) {
        for (auto iter=beg; iter<end; iter++) {
            _process_one_line(*iter, identifiedPeptides, outputFile, std::ref(m));
        }
        return;
    }
    auto mid = beg + len / 2;
    auto future = std::async(std::launch::async, parallel_match, mid, end, identifiedPeptides, outputFile, std::ref(m));
    parallel_match(beg, mid, identifiedPeptides, outputFile, std::ref(m));
    future.get();
}

int main() {
    // testing
    test_hamming1_search();
    std::mutex m;
    std::cout << "passed test!" << std::endl;
    char identified_filename[] = "/data/rqiao/wild_type_match/protein_origin.cc.txt";
    char output_filename[] = "/data/rqiao/wild_type_match/wildtype_matched.peptides.txt";
    char denovo_filename[] = "/data/rqiao/wild_type_match/SB_denovo.cc.txt";

    // read in identified peptides;
    auto identifiedPeptides = std::vector<std::vector<std::string>> ();

    readIdentifiedFile(identified_filename, identifiedPeptides);
    std::cout << identifiedPeptides.size() << " identified peptides" << std::endl;
    std::cout << identifiedPeptides[1].size() << std::endl;

    // match wild-type and write to output

    std::ifstream denovo_file(denovo_filename);
    std::string denovo_peptides, line_string;
    std::vector<std::string> denovo_peptide_list;
    std::vector<std::string> pattern;

    std::ofstream outputFile(output_filename);

    std::vector<std::future<void>> futures;

    while (getline(denovo_file, denovo_peptides)) {
        denovo_peptide_list.push_back(denovo_peptides);
    }

    parallel_match(denovo_peptide_list.begin(), denovo_peptide_list.end(), &identifiedPeptides, &outputFile, std::ref(m));

//    while (getline(denovo_file, denovo_peptides)) {
//        futures.push_back(std::async(std::launch::async, _process_one_line,
//                denovo_peptides, &identifiedPeptides, &outputFile, std::ref(m)));
//    }
//    std::cout << futures.size() << " peptides read" << std::endl;
//    for (auto &async_job : futures) {
//        async_job.get();
//    }
    outputFile.close();
    denovo_file.close();
    return 0;
}
