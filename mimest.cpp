/*
author: Nobutaka Kim
disclaimer: code uses snippets under 8 lines from Stackoverflow and other sources without citations
General algorithm from https://github.com/diogoff/unlabelled-event-logs
Use of OpenMP and Cuda to follow
*/
#include <cmath>
#include <vector>
#include <map>
#include <string>
#include <set>
#include <iostream>

void normalize(std::map<char, double> &gmNestedMap){
    double rowsum = 0.0;
    for(auto iter = gmNestedMap.begin(); iter!= gmNestedMap.end(); ++iter){
        char symbol = iter->first;
        double prob = iter->second;

        rowsum += prob; 
    }
    if(rowsum > 0.0){
        for(auto iter = gmNestedMap.begin(); iter!= gmNestedMap.end(); ++iter){
            gmNestedMap[iter->first] = (gmNestedMap[iter->first] / rowsum); 
        }
    }    
}


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
        std::map<char, int> M;
        //src sequence TBD
        std::vector<char> s;
        //separate source sequences (y^{(k)} in the paper)
        std::map<int, std::vector<char>> y;

        Model(std::vector<char> xIn, std::vector<char> DIn){
            x = xIn;
            N = x.size();
            //build D, set of symbols
            D.push_back(BEGIN);
            //eliminate duplicates and sort vector x
            //https://stackoverflow.com/a/1041939/9481613
            std::set<char> s;
            unsigned size = x.size();
            for(unsigned i = 0; i < size; ++i){
                s.insert(x[i]);
            }
            for(auto itr: s){
                D.push_back(itr);
            }
            D.push_back(END);

            for(auto &elementA : D){
                for(auto &elementB : D){
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
                std::map<char,double> &nestedMap = iter->second;
                normalize(nestedMap);
            }
        }

        void printModel(){

        }

        size_t estimate(){
            std::vector<std::vector<char>> prevsseqs;
            std::cout << "Initializing source sequence..." << std::endl;
        

        }

        size_t estsources(std::map<char, std::map<char, double>> T){
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

        size_t estparams(){
            char a;
            
            M.clear();
            for(int i=0; i < D.size(); i++){
                a = D[i];
                std::map<char, double> nestedMap;
                M[a] = nestedMap;
            }


        }
};


