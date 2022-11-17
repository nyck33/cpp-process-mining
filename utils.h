#ifndef UTILS_H
#define UTILS_H

#include <iostream>
#include <map>
#include <vector>
#include <algorithm>
#include <set>
#include <sstream>


static std::map<std::string, double> normalizeProbs(std::map<std::string, double> &probs){
    double rowsum = 0.0;
    for(auto & iter : probs){
        std::string str = iter.first;
        double prob = iter.second;

        rowsum += prob;
    }
    if(rowsum > 0.0){
        for(auto & iter : probs){
            //gmNestedMap[iter->first] = (gmNestedMap[iter->first] / rowsum);
            //gmNestedMap[iter->first] /= rowsum;
            //https://stackoverflow.com/a/45904048/9481613
            probs.at(iter.first) = (iter.second / rowsum);
        }
    }

    return probs;
}

static std::map<char, double> normalize(std::map<char, double> &gmNestedMap){
    //todo: need new instance of map to return for immutable
    double rowsum = 0.0;
    for(auto & iter : gmNestedMap){
        char symbol = iter.first;
        double prob = iter.second;

        rowsum += prob;
    }
    if(rowsum > 0.0){
        for(auto & iter : gmNestedMap){
            //gmNestedMap[iter->first] = (gmNestedMap[iter->first] / rowsum);
            //gmNestedMap[iter->first] /= rowsum;
            //https://stackoverflow.com/a/45904048/9481613
            gmNestedMap.at(iter.first) = (iter.second / rowsum);
        }
    }

    return gmNestedMap;
}

static std::string seq2str(std::vector<char> seq){
    std::string str = "";
    for(auto elem : seq){
        str += elem;
    }
    return str;
}

//1.  compare map values and sort map by value
bool cmp(std::pair<std::string, double>& a, std::pair<std::string, double>& b){
    return a.second > b.second;
}

void sortByValue(std::map<std::string, double>& M){
    std::vector<std::pair<std::string, double>> A;

    //copy key-value pair from map to vector of pairs
    for(auto& it : M){
        A.emplace_back(it);
    }

    //sort with comparator algo
    sort(A.begin(), A.end(), cmp);

    for(auto& it : A){
        std::cout << it.first << " " << it.second << std::endl;
    }

}


//open file, read lines and make vector of symbols
std::vector<char> openFileAndMakeVector(std::string inputFile){
    //todo: make x from file
    //std::string inputFile = "sequence.txt";
    std::vector<char> x;
    std::fstream newfile;
    newfile.open(inputFile, std::ios::in);

    if(!newfile.is_open()){
        perror("error open");
        exit(EXIT_FAILURE);
    }
    else if(newfile.is_open()){
        std::string tp;
        while(getline(newfile, tp)){
            //std::cout << tp << std::endl;
            if(isalpha(tp.at(0))){
                x.push_back(tp.at(0));
            }
            //std::cout << tp.size() << std::endl;

        }
        newfile.close();
    }

    return x;

}

static bool compareVectors(std::vector<int> a, std::vector<int> b){
    //check that vectors are same size
    if(a.size() != b.size()){
        return false;
    }
    
    std::sort(a.begin(), a.end());
    std::sort(b.begin(), b.end());
    return (a==b);
}

static bool compareSToPrevSeqs(std::vector<int> s, std::vector<std::vector<int>> prevseqs){
    bool vecsMatch = false;
    for(const auto& vec : prevseqs){
        vecsMatch = compareVectors(s, vec);
        if(vecsMatch){
            return true;
        }
    }
    //none matched
    return false;
}


//count unique items in vector
//https://stackoverflow.com/a/40851514/9481613
size_t count_items(const std::vector<int>& vec){
    std::set<int> counter;
    for(const int & x : vec){
        counter.insert(x);
    }
    size_t vecSize = counter.size();
    return vecSize;
}

//padding values for printing model
//https://stackoverflow.com/a/35451505/9481613
std::string pad_right(std::string const& str, size_t s)
{
    if ( str.size() < s )
        return str + std::string(s-str.size(), ' ');
    else
        return str;
}

std::string pad_left(std::string const& str, size_t s)
{
    if ( str.size() < s )
        return std::string(s-str.size(), ' ') + str;
    else
        return str;
}

std::string doubleToStringFixedPrecision(double num){
    std::ostringstream streamObj3;

    streamObj3 << std::fixed;
    streamObj3 << std::setprecision(2);
    streamObj3 << num;
    return streamObj3.str();
}

double roundValue(double val){
    val = std::ceil(val * 100.0) / 100.0;
    return val;
}

#endif