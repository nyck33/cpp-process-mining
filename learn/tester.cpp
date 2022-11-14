#include <iostream>
#include <map>
#include <vector>
#include <algorithm>

using std::cout; using std::cin;
using std::endl; using std::map;
using std::string; using std::vector;

int main(){
    /*
    map<string, string> veggy_map = {{"1", "Yam",},
                                     {"2", "Pumpkin",},
                                     {"3", "Ginger",},
                                     {"4", "Melon",},
                                     {"5", "Beetroot",},
                                     {"6", "Spinach",}};
    */
    map<string, double> dist = 
    {{"B", 0.04}, 
    {"A", 0.04},
    {"ACD", 0.17},
    {"GHEF", 0.086 },
    {"ACDF", 0.304}};

    cout << "Unsorted - " << endl;
    for (const auto & [key, value] : dist) {
        cout << key << " : " << value << endl;
    }

    cout << "Sorted - " << endl;

    map<string, string> veggy_map2;

    for (const auto & [key, value] : dist) {
        veggy_map2.emplace(value, key);
    }

    for (const auto & [key, value] : veggy_map2) {
        cout << value << " : " << key << endl;
    }

    return EXIT_SUCCESS;
}

