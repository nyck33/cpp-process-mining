#include <cmath>
#include <vector>
#include <map>
#include <string>

void normalize(){}

class Model{
    public:
        char BEGIN='o';
        char END = 'x';

        std::vector<char> x; 
        int N= 0;
        std::vector<char> D;

        std::map<int, std::string> gM;

        std::map<char, int> M;

        std::vector<char> s;

        std::map<char, char> y;

        




}