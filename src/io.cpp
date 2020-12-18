//
//  io.cpp
//

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include "io.hpp"
#include "basic_types.h"

std::vector<std::string> split(const std::string& s, char delimiter) {
    std::vector<std::string> tokens;
    std::string token;
    std::istringstream tokenStream(s);
    while (getline(tokenStream, token, delimiter)) {
        tokens.push_back(token);
    }
    return tokens;
}

int write_itemsets(std::string filename, VECTOR_TYPE** A_list,
                   VECTOR_TYPE** B_list, int list_size, int n, int m,
                   int *colOriginal, std::vector<std::string>attributes) {
    std::ofstream wf(filename);
    if (wf.is_open()) {
        for (int cid = 0; cid < list_size; cid++) {
            for (int i = 0; i < m; i++)
               if ((B_list[cid][i>>6]) & (VECTOR_1 << (i % 64)))
                    wf << attributes[colOriginal[i]] << " ";
            int s = 0;
            for (int i = 0; i < n; i++)
               if ((A_list[cid][i>>6]) & (VECTOR_1 << (i % 64)))
                    s++;
            wf << "; " << s << std::endl;
        }
        wf.close();
    } else {
        std::cout << "Cannot open file!" << std::endl;
        return 1;
    }
    return 0;
}

int read_itemsets(std::string filename) {
    std::vector<std::vector<std::string> > itemsets;
    std::vector<int> supports;
    std::string line;
    std::ifstream rf(filename, std::ios::in);
    if (!rf) {
       std::cout << "Cannot open file!" << std::endl;
       return 1;
    }
    while (getline(rf, line)) {
        std::vector<std::string> l = split(line, ';');
        assert(l.size() == 2);
        itemsets.push_back(split(l[0], ' '));
        supports.push_back(stoi(l[1]));
    }
    for (int i = 0; i <  itemsets.size(); i++) {
        std::vector<std::string> itemset = itemsets[i];
        for (int j = 0; j < itemset.size(); j++) {
            std::cout << itemset[j] << "  ";
        }
        std::cout << supports[i] << std::endl;
    }
    rf.close();
    return 0;
}

int write_summary(std::string filename, double* time_splits, int* level_size,
                  int* tree_size_array, int m) {
    std::ofstream wf(filename);
    if (wf.is_open()) {
        wf << "level;";
        for (int i = 0; i < m - 1; i++) {
            wf << i << ";";
        }
        wf << m - 1 << std::endl;
        wf << "n_itemsets;";
        for (int i = 0; i < m - 1; i++) {
            wf << level_size[i] << ";";
        }
        wf << level_size[m - 1] << std::endl;
        wf << "time;";
        for (int i = 0; i < m - 1; i++) {
            wf << time_splits[i] << ";";
        }
        wf << time_splits[m - 1] << std::endl;
        wf << "n_nodes;";
        for (int i = 0; i < m - 1; i++) {
            wf << tree_size_array[i] << ";";
        }
        wf << tree_size_array[m - 1] << std::endl;
        wf.close();
    } else {
        std::cout << "Cannot open file!" << std::endl;
        return 1;
    }
    return 1;
}




