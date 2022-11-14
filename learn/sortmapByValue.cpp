#include <iostream>
#include <map>
#include <vector>
#include <algorithm>


void sort(std::map<std::string, double> & M){

    //declare multimap
    std::multimap<double, std::string> MM;

    //insert each kv pair from M to MM as (vk) pairs
    for(auto& it: M){
        MM.insert({it.second, it.first});
    }

    //Print multimap
    for(auto& it : MM){
        std::cout << it.second << ' ' << it.first << std::endl;

    }
}

bool cmp(std::pair<std::string, double>& a,
std::pair<std::string, double>& b){
    return a.second > b.second;
}

void sortWithCmp(std::map<std::string, double>& M){
    std::vector<std::pair<std::string, double>> A;

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
int main(){
    //declare map
    std::map<std::string, double> dist = 
    {{"B", 0.04}, 
    {"A", 0.04},
    {"ACD", 0.17},
    {"GHEF", 0.086 },
    {"ACDF", 0.304}};

    sort(dist);

    std::cout << "\n" << std::endl;

    //use cmp sort
    sortWithCmp(dist);
    return 0;
}