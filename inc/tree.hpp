#include <stdlib.h>
#include <unordered_map>
#include <utility>

#ifndef INC_TREE_HPP_
#define INC_TREE_HPP_
struct node{
    bool isEndOfItemset;  // the end of the itemset
    std::unordered_map<int, node*> child;
};
std::pair<int, int>  checkTrie(node* cur_node);
node* insert(VECTOR_TYPE* itemset, node* &root, int m);
node* getNewNode();
void printTree();
size_t tree_size(node* cur_node);
size_t n_nodes(node* cur_node);
#endif  // INC_TREE_HPP_
