#ifndef UTILS_H
#define UTILS_H

#include <iostream>
#include <map>
#include <vector>
#include <algorithm>
#include <set>

//open file, read lines and make vector of symbols
std::vector<char> openFileAndMakeVector(std::string inputFile){
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
            std::cout << tp << std::endl;
            //std::cout << tp.size() << std::endl;

        }
        newfile.close();
    }

    return x;

}

//1.  compare map values and sort map by value
bool cmp(std::pair<char, double>& a, std::pair<char, double>& b){
    return a.second > b.second;
}

void sortByValue(std::map<char, double>& M){
    std::vector<std::pair<char, double>> A;

    //copy key-value pair from map to vector of pairs
    for(auto& it : M){
        A.push_back(it);
    }

    //sort with comparator algo
    sort(A.begin(), A.end(), cmp);

    for(auto& it : A){
        std::cout << it.first << " " << it.second << std::endl;
    }

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
    for(auto vec : prevseqs){
        vecsMatch = compareVectors(s, vec);
        if(vecsMatch == true){
            return true;
        }
    }
    return false;
}


static void normalize(std::map<char, double> &gmNestedMap){
    double rowsum = 0.0;
    for(auto iter = gmNestedMap.begin(); iter!= gmNestedMap.end(); ++iter){
        char symbol = iter->first;
        double prob = iter->second;

        rowsum += prob; 
    }
    if(rowsum > 0.0){
        for(auto iter = gmNestedMap.begin(); iter!= gmNestedMap.end(); ++iter){
            gmNestedMap[iter->first] = (gmNestedMap[iter->first] / rowsum); 
        }
    }    
}

static std::string seq2str(std::vector<char> seq){
    std::string str = "";
    for(auto elem : seq){
        str += elem;
    }
    return str;
}

//count unique items in vector
//https://stackoverflow.com/a/40851514/9481613
int count_items(std::vector<int> vec){
    std::set<int> counter;
    for(int & i : vec){
        counter.insert(i);
    }
    int vecSize = counter.size();
    return vecSize;
}


#endif