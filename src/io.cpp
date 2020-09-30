//
//  io.cpp
//  test
//
//  Created by Tatiana Makhalova on 01/07/2020.
//  Copyright © 2020 Tatiana Makhalova. All rights reserved.
//

#include "io.hpp"
#include "basic_types.h"
#include<iostream>
#include<fstream>
#include <sstream>
#include <string>
#include <vector>


using namespace std;

vector<string> split(const string& s, char delimiter)
{
   vector<string> tokens;
   string token;
   istringstream tokenStream(s);
   while (getline(tokenStream, token, delimiter))
   {
      tokens.push_back(token);
   }
   return tokens;
}


int write_itemsets(string filename, VECTOR_TYPE** A_list, VECTOR_TYPE** B_list, int list_size, int n, int m, int *colOriginal, vector<string>attributes){
    
    ofstream wf(filename);
    if(wf.is_open()) {
        for(int cid = 0; cid < list_size; cid++){
            for(int i = 0; i < m; i++)
                if((B_list[cid][i>>6])&(VECTOR_1<<(i%64)))
                    wf << attributes[colOriginal[i]] << " ";
            int s = 0;
            for(int i = 0; i < n; i++ )
                if((A_list[cid][i>>6])&(VECTOR_1<<(i%64)))
                    s++;
            wf << "; " << s << endl;
        }
        wf.close();
    }else{
        cout << "Cannot open file!" << endl;
        return 1;
    }
    return 0;
}

int read_itemsets(string filename){
    
    ifstream rf(filename, ios::in);
    if(!rf) {
        cout << "Cannot open file!" << endl;
        return 1;
    }
    vector<vector<string> > itemsets;
    vector<int> supports;
    string line;
    
    while(getline(rf,line)) {
        vector<string> l = split(line, ';');
        assert (l.size()==2);
        itemsets.push_back(split(l[0], ' '));
        supports.push_back(stoi(l[1]));
    }
    
//    output
    for (int i = 0; i <  itemsets.size(); i++){
        vector<string> itemset = itemsets[i];
        for(int j = 0; j < itemset.size(); j++){
            cout<<itemset[j] << "  ";
        }
        cout << supports[i] << endl;
    }
    rf.close();
   return 0;
}


int write_summary(string filename, double* time_splits, int* level_size, int* tree_size_array, int m){

    ofstream wf(filename);
    if(wf.is_open()){
        wf << "level;";
        for(int i = 0; i < m - 1; i++){
            wf << i << ";";
        }
        wf << m - 1 << endl;
        wf << "n_itemsets;";
        for(int i = 0; i < m - 1; i++){
            wf << level_size[i] << ";";
        }
        wf << level_size[m - 1] << endl;
        wf << "time;";
        for(int i = 0; i < m - 1; i++){
            wf << time_splits[i] << ";";
        }
        wf << time_splits[m - 1] << endl;
        wf << "n_nodes;";
        for(int i = 0; i < m - 1; i++){
            wf << tree_size_array[i] << ";";
        }
        wf << tree_size_array[m - 1] << endl;
        wf.close();
    }else{
        cout << "Cannot open file!" << endl;
        return 1;
    }
    return 1;
}




