//
// Created by rui on 11/8/18.
//

#ifndef WILD_TYPE_MATCH_STRING_SPLIT_H
#define WILD_TYPE_MATCH_STRING_SPLIT_H
#include <string>
#include <sstream>
#include <vector>
#include <iterator>

template<typename Out>
void split(const std::string &s, char delim, Out result) {
    std::stringstream ss(s);
    std::string item;
    while (std::getline(ss, item, delim)) {
        *(result++) = item;
    }
}

std::vector<std::string> split(const std::string &s, char delim) {
    std::vector<std::string> elems;
    split(s, delim, std::back_inserter(elems));
    return elems;
}

void join(const std::vector<std::string> &v, char c, std::string &s) {
    s.clear();
    for (auto p = v.begin(); p != v.end(); ++p) {
        s += *p;
        if (p != v.end() - 1)
            s += c;
    }
}

std::string join(const std::vector<std::string> &v, char c) {
    std::string s;
    s.clear();
    for (auto p = v.begin(); p != v.end(); ++p) {
        s += *p;
        if (p != v.end() - 1)
            s += c;
    }
    return s;
}

std::string join_subvector(const std::vector<std::string> &v, char c, unsigned long start_index, unsigned long length) {
    std::string s;
    s.clear();
    for (unsigned long i = start_index; i < start_index + length; i++) {
        s += v[i];
        if (i != start_index + length - 1)
            s += c;
    }
    return s;
}

#endif //WILD_TYPE_MATCH_STRING_SPLIT_H
