//
// Created by nycki_gq3buqc on 11/17/2022.
//

#ifndef CPP_PROCESS_MINING_MIMEST_H
#define CPP_PROCESS_MINING_MIMEST_H

#include <cmath>
#include <utility>
#include <vector>
#include <map>
#include <string>
#include <set>
#include <iostream>
#include <algorithm>
#include <fstream>
#include <iterator>
#include <cctype>
#include <iomanip>
#include <sstream>

struct sourcesRetStruct{
    std::vector<int> s;
    std::map<int, std::vector<char>> y;

};

struct estimateRetVal{
    size_t K;
    std::map<char, std::map<char, double>> M;

};



#endif //CPP_PROCESS_MINING_MIMEST_H
