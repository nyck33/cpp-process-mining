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

#include "utils.h"
#include "mimest.h"




int main(){
    char BEGIN='o';
    char END = 'x';
    char curr;
    //symbol sequence
    std::vector<char> x;

    std::string inputFile = "/home/nobu/Documents/ProcessMining/cpp-process-mining/mocksequence.txt";
    //std::string inputFile = "/home/nobu/Documents/ProcessMining/cpp-process-mining/learn/sequence.txt";
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
    gM = normalizeGM(intsD, gM);

    //transition matrix
    std::vector<std::vector<double>> M;
    //src sequence TBD
    std::vector<int> s;
    //separate source sequences (y^{(k)} in the paper)
    std::map<int, std::vector<char>> y;
    /////////////////////////////////////////////////
    //////call init funcs




    //estimate model
    estimateRetVal retVal = estimate(x,s,gM, M, D, y, N);
    size_t K = retVal.K;
    M = retVal.M;

    bool modelCorrect = checkmodel(x, y, s);
    std::cout << "model is correct: " << modelCorrect << std::endl;
    //print transition mat M
    printModel(M, D);

    std::map<std::string, double> probs = seqprobs(y);

    sortByValue(probs);

    std::cout << "Total number of sources: " << K << std::endl;

  
}


