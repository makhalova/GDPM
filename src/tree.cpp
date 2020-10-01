//
//  tree.cpp
//


#include <stdlib.h>
#include <iostream>
#include "basic_types.h"
#include "tree.hpp"

using namespace std;

node* getNewNode(){
    node* new_node = new node;
    new_node->isEndOfItemset= false;
    return new_node;
}

node * insert(VECTOR_TYPE *itemset, node* &root, int m){
    
    if (root == nullptr)
        root = getNewNode();
    
    int i = 0;
    bool last_exists = false;
    node *c_node = root;
    while (i < m){
        if (itemset[i>>6]&(VECTOR_1<<(i % 64))){ // if the item i is exist in 'itemset'
            if(c_node->child.find(i) == c_node->child.end()){
                c_node->child[i] = getNewNode();
                last_exists = false;
            }
            else{
                if (c_node->child[i]->isEndOfItemset)
                    last_exists = true;
                else
                    last_exists = false;
            }
            c_node = c_node->child[i];
        }
        i++;
    }
    if (last_exists)
        return nullptr;
    else
        c_node->isEndOfItemset = true;
    return c_node;
}

pair<int,int> checkTrie(node* cur_node){
    pair <int,int> n_nodes_init(0,0);
    if (cur_node->isEndOfItemset)
        n_nodes_init.first = 1; // true nodes
    else
        n_nodes_init.second = 1;

    for (pair<int, node*> element : cur_node->child){
        if (element.second != nullptr){
            pair<int,int>n_nodes = checkTrie(element.second);
            n_nodes_init.first += n_nodes.first;
            n_nodes_init.second += n_nodes.second;
        }
    }
    return n_nodes_init;
}


size_t tree_size(node* cur_node){
    if (cur_node->child.size() == 0)
        return sizeof(node);

    size_t total_size = 0;
    for (pair<int, node*> element : cur_node->child){
        total_size += tree_size(element.second);
    }
    return total_size;
}

size_t n_nodes(node* cur_node){
    if (cur_node->child.size() == 0)
        return 1;
    
    size_t total_size = 0;
    for (pair<int, node*> element : cur_node->child){
        total_size += n_nodes(element.second);
    }
    return total_size;
}
