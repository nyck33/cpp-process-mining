/*
https://www.geeksforgeeks.org/sorting-a-map-by-value-in-c-stl/
*/

#include <iostream>
#include <map>
#include <vector>
#include <algorithm>

//1.  compare map values and sort map by value
static bool cmp(std::pair<char, double>& a,
std::pair<char, double>& b){
    return a.second > b.second;
}

static void sortByValue(std::map<char, double>& M){
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

/*
//driver for sort map by value
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

