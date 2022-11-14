#include <cmath>
#include <vector>
#include <map>
#include <string>
#include <set>
#include <iostream>
#include <algorithm>
#include <fstream>

int main(){

    std::string inputFile = "sequence.txt";
    std::vector<char> x;
    std::fstream newfile;
    newfile.open(inputFile, std::ios::in);
    if(newfile.is_open()){
        std::string tp;
        while(getline(newfile, tp)){
            std::cout << tp << std::endl;
            std::cout << tp.size() << std::endl;

        }
        newfile.close();
    }

    
    
}