/*
author: Nobutaka Kim
disclaimer: code uses snippets under 8 lines from Stackoverflow and other sources without citations
General algorithm from https://github.com/diogoff/unlabelled-event-logs
Use of OpenMP and Cuda to follow

test with mimgen.py created sequence.text first then try porting
*/
#include <cmath>
#include <vector>
#include <map>
#include <string>
#include <set>
#include <iostream>
#include <algorithm>
#include <fstream>
#include <iterator>

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
	#   //char: {char : double}}
        std::map<char, std::map<char, double>> gM;
        //transition matrix
        std::map<char, std::map<char, double>> M;
        //src sequence TBD
        std::vector<int> s;
        //separate source sequences (y^{(k)} in the paper)
        std::map<int, std::vector<char>> y;

        Model(std::vector<char> xIn){
            x = xIn;
            N = x.size();
            //build D, set of symbols
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
            for(std::set<char>::iterator iter = sortedSetX.begin(); iter != sortedSetX.end(); iter++){
                D.push_back(*iter);
            }
            D.push_back(END);
            //init gM to 0.0
            for(auto elementA : D){
                for(auto elementB : D){
                    gM[elementA][elementB] = 0.0;
                }
            }
            //build gM
            int n;
            for(n=0; n < (N-1); ++n){
                char aChar = x[n];
                char bChar = x[n+1];
                gM[aChar][bChar] += 1.0;
            }
            for(auto iter = gM.begin(); iter!= gM.end(); ++iter){
                //make sure nestedMap is a reference 
                //std::map<char,double> &nestedMap = iter->second;
                //https://www.tutorialspoint.com/differences-between-pass-by-value-and-pass-by-reference-in-cplusplus
                normalize(iter->second);
            }
        }

        void printModel( std::map<char, std::map<char, double>> T){
            //https://stackoverflow.com/a/15553559/9481613
            std::cout << "Transition Matrix M:" << std::endl;
            for(auto iter = T.begin(); iter != T.end(); iter++){
                std::cout << iter->first << ": " << std::endl;
                std::map<char, double> nestedMap = iter->second;
                for(auto iter2 = nestedMap.begin(); iter2 != nestedMap.end(); iter2++){
                    std::cout << iter2->first << ": " << iter2->second << std::endl;
                }

            }
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
            char xn, b, bnext;

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
                            continue;
                        }
                    }
                    int vecsize = y[k].size();
                    char a = y[k][vecsize-1];
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
                for(int j = 0; j < D.size(); j++){
                    b = D[j];
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
            char a, b, k;
            
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
            for(auto a: D){
                normalize(M[a]);
            }
        }
        
        //computes the probability distribution for the different sequences produced by this model (p(z) or q(z) in the paper)
        std::map<char, double> seqprobs(){
            std::map<char, double> probs;
            for(auto & iter : y){
                std::string z = seq2str(iter.second);
                if(probs.find(z[0]) != probs.end()){
                    probs[z[0]] += 1.0;
                }else{
                    probs[z[0]] = 1.0;
                }

            }
            normalize(probs);
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

    std::string inputFile = "sequence.txt";
    std::fstream newfile;
    newfile.open(inputFile, std::ios::in);

    if(!newfile.is_open()){
        perror("error open");
        exit(EXIT_FAILURE);
    }
    else if(newfile.is_open()){
        std::string tp;
        //trim: https://stackoverflow.com/a/216883/9481613
        while(getline(newfile, tp)){
            std::cout << tp.at(0) << tp[0] << std::endl;
            curr = tp[0];
            x.push_back(curr);
        }
        newfile.close();
    }
    std::cout << "symbol sequence: " << seq2str(x) << std::endl;

    std::cout << x.size() << " symbols" << std::endl;

    //create model
    Model m = Model(x);

    //estimate model
    size_t K = m.estimate();

    bool modelCorrect = m.checkmodel();
    //print transition mat M
    m.printModel(m.M);

    std::map<char, double> probs = m.seqprobs(); 

    sortByValue(probs);

    std::cout << "Total number of sources: " << K << std::endl;

  
}


