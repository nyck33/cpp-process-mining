/*
https://www.geeksforgeeks.org/sorting-a-map-by-value-in-c-stl/
*/

#include <iostream>
#include <map>
#include <vector>
#include <algorithm>

bool cmp(std::pair<std::string, double>& a,
std::pair<std::string, double>& b){
    return a.second > b.second;
}

void sortByValue(std::map<std::string, double>& M){
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

/*
int main(){
    //declare map
    std::map<std::string, double> dist = 
    {{"B", 0.04}, 
    {"A", 0.04},
    {"ACD", 0.17},
    {"GHEF", 0.086 },
    {"ACDF", 0.304}};


    std::cout << "\n" << std::endl;

    //use cmp sort
    sortByValue(dist);
    return 0;
}
*/