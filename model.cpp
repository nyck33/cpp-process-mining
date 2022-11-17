//
// Created by nycki_gq3buqc on 11/16/2022.
//
//symbol sequence
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

std::string inputFile = "/home/nobu/Documents/ProcessMining/cpp-process-mining/mocksequence.txt";
//std::string inputFile = "/home/nobu/Documents/ProcessMining/cpp-process-mining/learn/sequence.txt";
//std::string inputFile = "sequence.txt";
