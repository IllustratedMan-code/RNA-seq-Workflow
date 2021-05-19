#include <iostream>
#include <fstream>
#include <string>
#include <vector>

int main(){
    std::ifstream in;
    in.open("GeneReference/genelist.fasta");
    std::string line;
    std::ofstream out;
    out.open("GeneReference/listofgenes.csv");
    if (in.is_open()){
        while(getline(in, line)){
            if (line[0] == '>'){
                line.erase(0, 1);
                out << line << std::endl;
            }
        }
    }
}
