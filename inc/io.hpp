//
//  io.hpp
//

#ifndef INC_IO_HPP_
#define INC_IO_HPP_
#include <stdio.h>
#include <vector>
#include <string>
#include <iostream>
#include "basic_types.h"

int write_binary_itemsets(std::string filename, VECTOR_TYPE** A_list,
                          VECTOR_TYPE** B_list, int list_size, int n, int m,
                          int* colOriginal);
int read_binary_itemsets(std::string filename);
int write_itemsets(std::string filename, VECTOR_TYPE** A_list,
                   VECTOR_TYPE** B_list, int list_size, int n, int m,
                   int* colOriginal, std::vector<std::string>attributes);
int read_itemsets(std::string filename);
int write_summary(std::string filename, double* time_splits, int* level_size,
                  int* tree_size_array, int m);

#endif   // INC_IO_HPP_

