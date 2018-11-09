//
// Created by rui on 11/8/18.
//

#ifndef WILD_TYPE_MATCH_HAMMING1_SEARCH_H
#define WILD_TYPE_MATCH_HAMMING1_SEARCH_H

#include <vector>
#include <string>
#include <stdexcept>
#include <map>

struct matchLocations {
    bool exactMatch;
    std::vector<unsigned long> startIndex;
};


matchLocations hamming1_search(const std::vector<std::string> &pattern, const std::vector<std::string> &text) {
    // if there is an exact match, stop search
    unsigned long m = pattern.size();
    int R0=0, R1=0, S=0, shR0;
    int mask = 1 << (m - 1);
    if (m > 30) {
        throw std::runtime_error("pattern length should be less than 31");
    }

    auto result = matchLocations {
        false,
        {},
    };

    std::map<std::string, int> sMap;
    int index = 0;
    // initialize sMap
    for (const std::string &c : pattern) {
        if (sMap.find(c) == sMap.end()) {
            sMap.insert(std::make_pair(c, 1 << index));
        } else {
            sMap[c] |= 1 << index;
        }
        index++;
    }

    index = 0;
    for (const std::string &c : text) {
        if (sMap.find(c) == sMap.end()) {
            S = 0;
        } else {
            S = sMap[c];
        }
        shR0 = (R0 << 1) | 1;
        R0 = shR0 & S;
        R1 = ( ((R1 << 1) | 1) & S ) | shR0;
        if (R0 & mask) {
            result.exactMatch = true;
            // found an exact match, immediately return.
            return result;
        } else if (R1 & mask) {
            result.startIndex.push_back(index - m + 1);
        }
        index++;
    };

    return result;
}

#endif //WILD_TYPE_MATCH_HAMMING1_SEARCH_H
