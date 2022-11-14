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
#include <algorithm>

static bool compareSToPrevSeqs(std::vector<int> s, std::vector<std::vector<int>> prevseqs){
    bool vecsMatch = false;
    for(auto vec : prevseqs){
        vecsMatch = compareVectors(s, vec);
        if(vecsMatch == true){
            return true;
        }
    }
    return false;
}

static bool compareVectors(std::vector<int> a, std::vector<int> b){
    //check that vectors are same size
    if(a.size() != b.size()){
        return false;
    }
    
    std::sort(a.begin(), a.end());
    std::sort(b.begin(), b.end());
    return (a==b);
}

static void normalize(std::map<char, double> &gmNestedMap){
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

static std::string seq2str(std::vector<char> seq){
    std::string str = "";
    for(auto elem : seq){
        str += elem;
    }
    return str;
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
        std::map<char, std::map<char, double>> M;
        //src sequence TBD
        std::vector<int> s;
        //separate source sequences (y^{(k)} in the paper)
        std::map<int, std::vector<char>> y;

        Model(std::vector<char> xIn, std::vector<char> DIn){
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
            //sort
            sort(sortedSetX.begin(), sortedSetX.end());
            //push into vector
            for(auto itr: s){
                D.push_back(itr);
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

        void printModel(){

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

            return std::set(s.begin(), s.end()).size();    

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
            char a, b, k;
            
            M.clear();
            for(int i=0; i < D.size(); i++){
                a = D[i];
                std::map<char, double> nestedMap;
                for(int j=0; j < D.size(); j++){
                    b = D[j];
                    M[a][b] = 0.0; 
                }
            }
            for(auto iter = y.begin(); iter!= y.end(); ++iter){
                k = iter->first;
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
        std::map<char, double> computePDistForSequences(){
            std::map<char, double> probs;
            for(auto iter = y.begin(); iter != y.end(); ++iter){
                
            }



            return probs;

        }

        //checks that it is possible to recover the symbol sequence x from the separate sequences y (sanity check)
        bool checkmodel(){
            int sn;
            char xn;

            std::vector<char> x2;
            std::map<int, int> pos;
            for(auto iter = y.begin(); iter!= y.end(); ++iter){
                pos[iter->first] = -1;
            }
            for(int n = 0; n < s.size(); n++){
                sn = s[n];
                pos[sn] += 1;
                xn = y[sn][pos[sn]];
                x2.push_back(xn);
            }
            return (x2 == x);


        }


};


