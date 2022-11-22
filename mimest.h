//
// Created by nycki_gq3buqc on 11/17/2022.
//

#ifndef CPP_PROCESS_MINING_MIMEST_H
#define CPP_PROCESS_MINING_MIMEST_H

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

struct sourcesRetStruct{
    std::vector<int> s;
    std::map<int, std::vector<int>> y;


};

struct estimateRetVal{
    size_t K{};
    std::vector<std::vector<double>> M;
    std::map<int, std::vector<int>> intsY;


};

struct DdictStruct{
    std::map<char, int> Ddict;
    std::map<int, char> revDdict;

};

std::vector<int> translateXToInts(const std::vector<char>& x, std::map<char, int>& dDict){
    std::vector<int> intsX(x.size(), 0);
    //std::vector<int>::iterator it;
    #pragma omp parallel for default(none) shared(x, dDict, intsX)
    for(int i =0; i < x.size(); i++){
        intsX[i] = dDict[x[i]];
    }
    /*
    for(auto it: x){
        intsX.push_back(dDict[it]);
    }
    */
    return intsX;
}

DdictStruct makeDdict(std::vector<char> D){
    std::map<char, int> Ddict;
    std::map<int, char> revDdict;

    #pragma omp parallel for default(none) shared(revDdict, Ddict, D)
    for(int i = 0; i < D.size(); i++ ){
        Ddict[D[i]] = i;
        revDdict[i] = D[i];
    }
    DdictStruct dictStruct;
    dictStruct.Ddict = Ddict;
    dictStruct.revDdict = revDdict;

    return dictStruct;
}


std::vector<char> initializeD(std::vector<char> x, char begin, char end) {
    std::vector<char> D;
    D.push_back(begin);
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
    D.push_back(end);

    return D;
}

std::vector<std::vector<double>> initializeGM(std::vector<int>& D) {

    std::vector<double> innerVec(D.size(), 0.0);
    std::vector<std::vector<double>> gM(D.size(),innerVec);


    return gM;
}

std::vector<std::vector<double>> buildGM(std::vector<int> x,
                                         std::vector<std::vector<double>> gM,
                                               size_t N) {

    //int a, b;
    //std::map<char, std::map<char, double>> newGM;
    #pragma omp parallel for default(none) shared(N, x, gM)
    for (int n = 0; n < N - 1; n++) {
        int a = x.at(n);
        int b = x[n + 1];
        gM.at(a).at(b) += 1.0;
    }

    return gM;
}

std::vector<std::vector<double>> normalizeGM(std::vector<std::vector<double>> &gM) {
    double rowsum;
    std::vector<double> rowsumsVec(gM.size(), 0);

    #pragma omp parallel for default(none) shared(gM, rowsumsVec)
    for(int i=0; i<gM.size(); i++)
    for (int j = 0; j<gM[i].size(); j++) {
        rowsumsVec[i] += gM[i][j];

        if(rowsumsVec[i] > 0) {
            for (int k = 0; k < gM[i].size(); k++) {
                gM[i][k] /= rowsumsVec[i];
            }
        }

    }

    return gM;

}

void printModel(std::vector<std::vector<double>> T, std::vector<char> D,
                std::map<char, int> Ddict) {
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
    //very left row names
    for (auto a: D) {
        outputStr = std::string(1, a);
        outputStr = pad_right(outputStr, padSize);
        std::cout << outputStr;
        //columns
        for (auto b: D) {
            int aIdx, bIdx;
            aIdx = Ddict[a];
            bIdx = Ddict[b];
            if (T.at(aIdx).at(bIdx) == 0.0) {
                outputStr = pad_right(blank, padSize);
                std::cout << outputStr;
            } else {
                val = T.at(aIdx).at(bIdx);
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
std::map<std::string, double> seqprobs(std::map<int, std::vector<char>> &y){
    std::map<std::string, double> probs;

    #pragma omp parallel for default(none) shared(probs, y)
    for(int i=0; i<y.size(); i++){
        if(probs.find(seq2str(y[i]))!= probs.end()){
            probs[seq2str(y[i])] += 1.0;
        }else{
            probs[seq2str(y[i])] = 1.0;
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

sourcesRetStruct estsources(std::vector<int> x,
                            const std::vector<int>& D,
                            size_t N,
                            std::vector<std::vector<double>> T, int BEGIN, int END) {
    //double pmax, p, pnext;
    //int sn;
    // char xn, b, bnext, a;

    sourcesRetStruct retVals;

    std::vector<int> s;
    std::map<int, std::vector<int>> y;
    std::set<int> active;
    bool xnInYk = false;

    for (int n = 0; n < N; n++) {
        int xn = x[n];
        double pmax = 0.0;
        int sn = -1;
        for (int k: active) {
            for (auto elem: y[k]) {
                if (xn == elem) {
                    //if xn in y[k]
                    xnInYk = true;
                    break;
                }
            }
            if(xnInYk){
                xnInYk = false;
                continue;
            }
            int vecsize = y[k].size();
            // std::map<int, std::vector<char>> y;
            int a = y[k][vecsize - 1];
            int b = xn;
            double p = T[a][b];
            if (p > pmax) {
                sn = k;
                pmax = p;
            }
        }
        double TBeginXn = T[BEGIN][xn];
        if ((sn == -1) || (TBeginXn > pmax)) {
            sn = y.size() + 1;
            active.insert(sn);
            if(!y[sn].empty()){
                y[sn].clear();
            }
        }
        s.push_back(sn);
        y[sn].push_back(xn);
        double pnext = 0.0;
        int bnext = BEGIN;
        for (int b: D) {
            double Txnb = T[xn][b];
            if (Txnb > pnext) {
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
std::vector<std::vector<double>> estparams(
        std::vector<int>& D,
        std::map<int, std::vector<int>>& y, int BEGIN, int END) {
    //std::map<char, std::map<char, double>> M;
    //int k;
    std::vector<double>MRow(D.size(), 0.0);
    std::vector<std::vector<double>> M(D.size(),MRow);

    //int k = 0;

    //#pragma omp parallel for default(none) shared(y, BEGIN, M, END)
    for(int k = 1; k< y.size(); k++){
        int a = BEGIN;
        int temp = k;
        int b = y[k][0];

        M.at(a).at(b) += 1.0;

        for(int r=0; r < y[k].size()-1; r++){
            a = y[k][r];
            b = y[k][r+1];
            M[a][b] += 1.0;
        }
        a = y[k][y[k].size()-1];
        b = END;
        M[a][b] += 1.0;
    }

    M = normalizeGM(M);

    return M;
}


//estimate(x,s,gM, M, D, N);
estimateRetVal estimate(std::vector<int>& x,
                        std::vector<int> s,
                        std::vector<std::vector<double>>& gM,
                        std::vector<std::vector<double>> M,
                        std::vector<int>& D,
                        std::map<int, std::vector<int>> y,
                        size_t N, int BEGIN, int END, std::map<char, int> &dDict,
                        std::vector<char> &charD) {
    std::vector<std::vector<int>> prevsseqs;
    std::cout << "Initializing source sequence..." << std::endl;
    sourcesRetStruct srcsSY = estsources(x, D, N, gM, BEGIN, END);
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
        M = estparams(D,y, BEGIN, END);
        prevsseqs.push_back(s);
        std::cout << "# " << its << ": Computing source sequence..." << std::endl;
        srcsSY = estsources(x, D, N, M, BEGIN, END);
        y = srcsSY.y;
        s = srcsSY.s;

        printModel(M, charD, dDict);

    }
    //return num unique items
    estimateRetVal retValStruct;
    retValStruct.K = count_items(s);
    retValStruct.M = M;
    retValStruct.intsY = y;
    return retValStruct;

}



#endif //CPP_PROCESS_MINING_MIMEST_H
