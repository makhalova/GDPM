#include <stdlib.h> /* malloc, free, rand */
#include <unordered_map>



#ifndef tree_hpp
#define tree_hpp

struct node{
    
    bool isEndOfItemset; // the end of the itemset
    std::unordered_map<int, node*> child;
    
};

std::pair<int, int>  checkTrie(node* cur_node);

node * insert(VECTOR_TYPE* itemset, node* &root, int m); //node*& root,
node * getNewNode();

void printTree();

size_t tree_size(node* cur_node);
size_t n_nodes(node* cur_node);
    
#endif
