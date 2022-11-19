/*
author: Nobutaka Kim
disclaimer: code uses snippets under 8 lines from Stackoverflow and other sources without citations
General algorithm from https://github.com/diogoff/unlabelled-event-logs
Use of OpenMP and Cuda to follow

test with mimgen.py created sequence.text first then try porting
*/
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
#include <omp.h>

#include "utils.h"
#include "mimest.h"




int main(){
    char BEGIN='o';
    char END = 'x';
    char curr;
    //symbol sequence
    std::vector<char> x;

    //std::string inputFile = "/home/nobu/Documents/ProcessMining/cpp-process-mining/mocksequence.txt";
    std::string inputFile = "/home/nobu/Documents/ProcessMining/cpp-process-mining/learn/sequence.txt";
    x = openFileAndMakeVector(inputFile);

    std::cout << "symbol sequence: " << seq2str(x) << std::endl;

    std::cout << x.size() << " symbols" << std::endl;

    size_t N = x.size();

    //set of symbols in x
    std::vector<char> D;
    D = initializeD(x, BEGIN, END);
    DdictStruct dictStruct = makeDdict(D);
    std::map<char, int> dDict = dictStruct.Ddict;
    std::map<int, char> revDDict = dictStruct.revDdict;

    std::vector<int> intsX = translateXToInts(x, dDict);
    std::vector<int> intsD = translateXToInts(D, dDict);

    std::vector<std::vector<double>> gM = initializeGM(intsD);

    gM = buildGM(intsX, gM, N);
    gM = normalizeGM(gM);

    size_t numRows = gM.size();
    size_t numCols = gM[0].size();

    //transition matrix
    std::vector<std::vector<double>> M = initializeGM(intsD);
    //src sequence TBD
    std::vector<int> s;
    //separate source sequences (y^{(k)} in the paper)
    //each nested vec can differ in size
    std::map<int, std::vector<int>> y;

    //redefine begin and end with ints
    int BEGINInt = dDict[BEGIN];
    int ENDInt = dDict[END];
    /////////////////////////////////////////////////
    //////call init funcs

    //estimate model
    estimateRetVal retVal = estimate(intsX,s,gM, M, intsD, y, N, BEGINInt, ENDInt);
    size_t K = retVal.K;
    M = retVal.M;
    y = retVal.intsY;

    //translate y back to chars from int
    std::map<int, std::vector<char>> charY;

    #pragma omp parallel for
    for(int i=0; i<y.size(); i++){
        std::vector<char> subarr;
        for(int j=0; j < y[i].size(); j++){
            subarr.push_back(revDDict[y[i][j]]);
        }
        charY[i] = subarr;
    }


    //bool modelCorrect = checkmodel(intsX, y, s);
    //std::cout << "model is correct: " << modelCorrect << std::endl;
    //print transition mat M
    printModel(M, D, dDict);

    std::map<std::string, double> probs = seqprobs(charY);

    sortByValue(probs);

    std::cout << "Total number of sources: " << K << std::endl;

  
}


