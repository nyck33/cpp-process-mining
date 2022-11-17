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

char BEGIN='o';
char END = 'x';

std::vector<char> initializeD(std::vector<char> x) {
    std::vector<char> D;
    D.push_back(BEGIN);
    //eliminate duplicates and sort vector x
    //https://stackoverflow.com/a/1041939/9481613
    std::set<char> sortedSetX;
    unsigned size = x.size();
    //eliminate duplicates
    for (unsigned i = 0; i < size; ++i) {
        sortedSetX.insert(x[i]);
    }
    //sort: https://linuxhint.com/sorting-elements-cpp-set/
    for (char iter: sortedSetX) {
        if (std::isalpha(iter)) {
            D.push_back(iter);
        }
    }
    D.push_back(END);

    return D;
}


std::map<char, std::map<char, double>> initializeGM(const std::vector<char>& D) {
    std::map<char, std::map<char, double>> gM;
    //todo: init gM to 0.0
    for (auto elementA: D) {
        for (auto elementB: D) {
            gM[elementA][elementB] = 0.0;
        }
    }
    return gM;
}

std::map<char, std::map<char, double>> buildGM(std::vector<char> x,
                                               std::map<char, std::map<char, double>> gM,
                                               size_t N) {

    char a, b;
    //std::map<char, std::map<char, double>> newGM;

    for (int n = 0; n < N - 1; n++) {
        a = x[n];
        b = x[n + 1];
        gM[a][b] += 1.0;
    }

    return gM;
}

std::map<char, std::map<char, double>> normalizeGM(std::vector<char> D,
                                                   std::map<char, std::map<char, double>> gM) {
    std::map<char, std::map<char, double>> newGM;

    for (auto a: D) {
        newGM[a] = normalize(gM[a]);

    }

    return newGM;

}

void printModel(std::map<char, std::map<char, double>> T, std::vector<char> D) {
    double val;
    std::string outputStr, blank = "-";

    size_t padSize = 5;
    //https://stackoverflow.com/a/15553559/9481613
    std::cout << "Transition Matrix M:" << std::endl;
    //5 spaces
    outputStr = pad_right(blank, padSize);
    std::cout << outputStr;
    //top columns
    for (auto a: D) {
        outputStr = std::string(1, a);
        outputStr = pad_right(outputStr, padSize);
        std::cout << outputStr;
    }
    std::cout << std::endl;
    for (auto a: D) {
        outputStr = std::string(1, a);
        outputStr = pad_right(outputStr, padSize);
        std::cout << outputStr;

        for (auto b: D) {
            if (T[a][b] == 0.0) {
                outputStr = pad_right(blank, padSize);
                std::cout << outputStr;
            } else {
                val = T[a][b];
                outputStr = doubleToStringFixedPrecision(val);
                outputStr = pad_right(outputStr, padSize);
                std::cout << outputStr;
            }

        }
        std::cout << std::endl;
    }
    /*
    for(auto & iter : T){
        std::cout << iter.first << ": " << std::endl;
        std::map<char, double> nestedMap = iter.second;
        for(auto & iter2 : nestedMap){
            std::cout << iter2.first << ": " << iter2.second << std::endl;
        }
    }
    */
}

//computes the probability distribution for the different sequences produced by this model (p(z) or q(z) in the paper)
//        std::map<int, std::vector<char>> y;
std::map<std::string, double> seqprobs(const std::map<int, std::vector<char>>& y){
    std::map<std::string, double> probs;
    for(auto & iter : y){
        std::string z = seq2str(iter.second);
        if(probs.find(z) != probs.end()){
            probs[z] += 1.0;
        }else{
            probs[z] = 1.0;
        }

    }
    normalizeProbs(probs);
    return probs;
}

//checks that it is possible to recover the symbol sequence x from the separate sequences y (sanity check)
bool checkmodel(const std::vector<char>& x,
                std::map<int, std::vector<char>> y,
                std::vector<int> s){
    int sn;
    char xn;

    std::vector<char> x2;
    std::map<unsigned long, int> pos;
    for(auto & iter : y){
        pos[iter.first] = -1;
    }
    for(int n : s){
        sn = n;
        pos[sn] += 1;
        xn = y[sn][pos[sn]];
        x2.push_back(xn);
    }
    return (x2 == x);
}

sourcesRetStruct estsources(std::vector<char> x,
                            std::vector<char> D,
                            size_t N,
                            std::map<char, std::map<char, double>> T) {
    //double pmax, p, pnext;
    //int sn;
    // char xn, b, bnext, a;

    sourcesRetStruct retVals;

    std::vector<int> s;
    std::map<int, std::vector<char>> y;
    std::set<int> active;
    for (int n = 0; n < N; n++) {
        char xn = x[n];
        double pmax = 0.0;
        long sn = -1;
        for (int k: active) {
            for (auto elem: y[k]) {
                if (xn == elem) {
                    //if xn in y[k]
                    continue;
                }
            }
            int vecsize = y[k].size();
            // std::map<int, std::vector<char>> y;
            char a = y[k][vecsize - 1];
            char b = xn;
            double p = T[a][b];
            if (p > pmax) {
                sn = k;
                pmax = p;
            }
        }
        if ((sn == -1) || (T[BEGIN][xn] > pmax)) {
            sn = y.size() + 1;
            active.insert(sn);
            y[sn].clear();
        }
        s.push_back(sn);
        y[sn].push_back(xn);
        double pnext = 0.0;
        char bnext = BEGIN;
        for (char b: D) {
            if (T[xn][b] > pnext) {
                pnext = T[xn][b];
                bnext = b;
            }
        }
        if (bnext == END) {
            active.erase(sn);
        }
    }
    retVals.s = s;
    retVals.y = y;

    return retVals;
}

//return M
std::map<char, std::map<char, double>> estparams(
        const std::vector<char>& D,
        std::map<int, std::vector<char>> y) {
    //std::map<char, std::map<char, double>> M;
    //int k;

    std::map<char, std::map<char, double>> M;
    for (char a: D) {
        //std::map<char, double> nestedMap;
        for (char b: D) {

            M[a][b] = 0.0;
        }
    }
    for (auto &iter: y) {
        //std::map<int, std::vector<char>> y;
        int k = iter.first;
        char a = BEGIN;
        char b = y[k][0];
        M[a][b] += 1.0;

        for (int r = 0; r < (y[k].size() - 1); r++) {
            a = y[k][r];
            b = y[k][r + 1];
            M[a][b] += 1.0;
        }
        a = y[k][y.size() - 1];
        b = END;
        M[a][b] += 1.0;
    }
    for (auto aChar: D) {
        M[aChar] = normalize(M[aChar]);
    }
    return M;
}
//estimate(x,s,gM, M, D, N);
estimateRetVal estimate(const std::vector<char>& x,
                        std::vector<int> s,
                        std::map<char, std::map<char, double>> gM,
                        std::map<char, std::map<char, double>> M,
                        const std::vector<char>& D,
                        std::map<int, std::vector<char>> y,
                        size_t N) {
    std::vector<std::vector<int>> prevsseqs;
    std::cout << "Initializing source sequence..." << std::endl;
    sourcesRetStruct srcsSY = estsources(x, D, N, std::move(gM));
    s = srcsSY.s;
    y = srcsSY.y;
    int its = 0;
    //compare s to each vector in vector of vectors prevsseqs
    bool done = false;
    while (true) {
        done = compareSToPrevSeqs(s, prevsseqs);
        if(done){
            break;
        }
        its += 1;
        std::cout << "# " << its << ": Estimating parameters..." << std::endl;
        M = estparams(D,y);
        prevsseqs.push_back(s);
        std::cout << "# " << its << ": Computing source sequence..." << std::endl;
        srcsSY = estsources(x, D, N, M);
        y = srcsSY.y;
        s = srcsSY.s;

    }
    //return num unique items
    estimateRetVal retValStruct;
    retValStruct.K = count_items(s);
    retValStruct.M = M;
    return retValStruct;

}

int main(){
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
    //global model used to initialize M, (M^{+} in the paper)
    //# nested {state: {nextstate: prob, nextstate2: prob2,...}}
    //char: {char : double}}
    std::map<char, std::map<char, double>> gM;
    //transition matrix
    std::map<char, std::map<char, double>> M;
    //src sequence TBD
    std::vector<int> s;
    //separate source sequences (y^{(k)} in the paper)
    std::map<int, std::vector<char>> y;

    /////////////////////////////////////////////////
    //////call init funcs
    D = initializeD(x);
    gM = initializeGM(D);
    gM = buildGM(x, gM, N);
    gM = normalizeGM(D, gM);

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


