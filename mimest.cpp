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

class Model{
    public:
        char BEGIN='o';
        char END = 'x';
        //symbol sequence
        std::vector<char> x; 
        //len of x
        int N= 0;
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

        Model(std::vector<char> xIn){
            x = std::move(xIn);
            N = x.size();
            //todo: build D, set of symbols
            D.push_back(BEGIN);
            //eliminate duplicates and sort vector x
            //https://stackoverflow.com/a/1041939/9481613
            std::set<char> sortedSetX;
            unsigned size = x.size();
            //eliminate duplicates
            for(unsigned i = 0; i < size; ++i){
                sortedSetX.insert(x[i]);
            }
            //sort: https://linuxhint.com/sorting-elements-cpp-set/
            for(char iter : sortedSetX){
                if(std::isalpha(iter)){
                    D.push_back(iter);
                }
            }
            D.push_back(END);
            //todo: init gM to 0.0
            for(auto elementA : D){
                for(auto elementB : D){
                    gM[elementA][elementB] = 0.0;
                }
            }
            //todo: build gM, each time symbol appears in seq x, +1.0
            int n;
            for(n=0; n < (N-1); ++n){
                char aChar = x[n];
                char bChar = x[n+1];
                gM[aChar][bChar] += 1.0;
            }
            //todo: normailize nested map of gM
            for(auto & iter : gM){
                //make sure nestedMap is a reference 
                //std::map<char,double> &nestedMap = iter->second;
                //https://www.tutorialspoint.com/differences-between-pass-by-value-and-pass-by-reference-in-cplusplus
                normalize(iter.second);
            }
        }

        void printModel( std::map<char, std::map<char, double>> T){
            double val, roundedVal;
            std::string outputStr, blank = "-";

            size_t padSize =5;
            //https://stackoverflow.com/a/15553559/9481613
            std::cout << "Transition Matrix M:" << std::endl;
            //5 spaces
            outputStr = pad_right(blank, padSize);
            std::cout << outputStr;
            //top columns
            for(auto a : D){
                outputStr = std::string(1, a);
                outputStr = pad_right(outputStr, padSize);
                std::cout << outputStr;
            }
            std::cout << std::endl;
            for(auto a: D){
                outputStr = std::string(1,a);
                outputStr = pad_right(outputStr, padSize);
                std::cout << outputStr;

                for(auto b: D){
                    if(T[a][b] == 0.0){
                        outputStr = pad_right(blank, padSize);
                        std::cout << outputStr;
                    }else{
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

        size_t estimate(){
            std::vector<std::vector<int>> prevsseqs;
            std::cout << "Initializing source sequence..." << std::endl;
            estsources(gM);
            int its = 0;
            //compare s to each vector in vector of vectors prevsseqs
            while(!compareSToPrevSeqs(s,prevsseqs)){
                its += 1;
                std::cout << "# " << its << ": Estimating parameters..." << std::endl;
                estparams();
                prevsseqs.push_back(s);
                std::cout << "# " << its << ": Computing source sequence..." << std::endl;
                estsources(M);

            }   
            //return num unique items
            return count_items(s);

        }

        void estsources(std::map<char, std::map<char, double>> T){
            double pmax, p, pnext;
            int sn;
            char xn, b, bnext, a;

            s.clear();
            y.clear();
            std::set<int> active;
            for(int n=0; n<N; n++){
                xn = x[n];
                pmax = 0.0;
                sn = -1;
                for(auto k: active){
                    for(auto elem : y[k]){
                        if(xn == elem){
                            //if xn in y[k]
                            continue;
                        }
                    }
                    int vecsize = y[k].size();
                    a = y[k][vecsize-1];
                    b = xn;
                    p = T[a][b];
                    if(p > pmax){
                        sn = k;
                        pmax = p;
                    }
                }
                if((sn == -1) || (T[BEGIN][xn] > pmax)){
                    sn = y.size() + 1;
                    active.insert(sn);
                    y[sn].clear();
                }
                s.push_back(sn);
                y[sn].push_back(xn);
                pnext = 0.0;
                bnext = BEGIN;
                for(char j : D){
                    b = j;
                    if(T[xn][b] > pnext){
                        pnext = T[xn][b];
                        bnext = b;
                    }
                    if(bnext == END){
                        active.erase(sn);
                    }
                }   
            }
        }   

        void estparams(){
            char a, b;
            int k;
            
            M.clear();
            for(char i : D){
                a = i;
                std::map<char, double> nestedMap;
                for(char j : D){
                    b = j;
                    M[a][b] = 0.0; 
                }
            }
            for(auto & iter : y){
                //std::map<int, std::vector<char>> y;
                k = iter.first;
                a = BEGIN;
                b = y[k][0];

                for(int r=0; r < (y[k].size() - 1); r++){
                    a = y[k][r];
                    b = y[k][r+1];
                    M[a][b] += 1.0;
                }
                a = y[k][y.size()-1];
                b = END;
                M[a][b] += 1.0;
            }
            for(auto aChar: D){
                normalize(M[aChar]);
            }
        }
        
        //computes the probability distribution for the different sequences produced by this model (p(z) or q(z) in the paper)
        //        std::map<int, std::vector<char>> y;
        std::map<std::string, double> seqprobs(){
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
        bool checkmodel(){
            int sn;
            char xn;

            std::vector<char> x2;
            std::map<int, int> pos;
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
};

int main(){
    std::vector<char> x;
    char curr;

    //std::string inputFile = "/home/nobu/Documents/ProcessMining/cpp-process-mining/mocksequence.txt";
    std::string inputFile = "/home/nobu/Documents/ProcessMining/cpp-process-mining/learn/sequence.txt";
    //std::string inputFile = "sequence.txt";
    x = openFileAndMakeVector(inputFile);
    std::cout << "symbol sequence: " << seq2str(x) << std::endl;

    std::cout << x.size() << " symbols" << std::endl;

    //create model
    Model m = Model(x);

    //estimate model
    size_t K = m.estimate();

    bool modelCorrect = m.checkmodel();
    std::cout << "model is correct: " << modelCorrect << std::endl;
    //print transition mat M
    m.printModel(m.M);

    std::map<std::string, double> probs = m.seqprobs();

    sortByValue(probs);

    std::cout << "Total number of sources: " << K << std::endl;

  
}


