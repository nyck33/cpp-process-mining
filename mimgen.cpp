#include <cmath>
#include <vector>
#include <map>
#include <string>
#include <fstream>
#include <iostream>
#include <unistd.h>
#include <algorithm>

//func prototypes
void printseqs();

//global variables

char* getCmdOptions(char** begin, char** end, const std::string &option){
    char** itr = std::find(begin, end, option);
    if(itr != end && ++itr != end){
        return *itr;
    }
    return 0;
}

bool cmdOptionExists(char** begin, char** end, const std::string &option){
    //single word option ie. -h for help
    return std::find(begin, end, option) != end;
}

int main(int argc, char* argv[]){
    int ninstances, overlap;
    char* inputDistributionFile, outputSymbolFile;
    //read cmd line arguments: https://stackoverflow.com/a/868894/9481613
    //mimgen 20 5 example.txt sequence.txt
    //prog --instances 20 --max-overlapping 5 --inputDist example --outSymbol sequence 
    //prog -n 20 -m 5 -i example.txt -o sequence.txt

    char* numInstancesStr = getCmdOptions(argv, argv+argc, "-n");
    
    if(numInstancesStr){
        //https://stackoverflow.com/a/23567980/9481613
        ninstances = std::stoi(numInstancesStr);        
    }

    char* overlapStr = getCmdOptions(argv, argv+argc, "-m");

    if(overlapStr){
        overlap = std::stoi(overlapStr);
    }

}





        
        




