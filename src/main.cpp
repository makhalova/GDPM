//  GDPM (version 1.0)
//
//  Created by Tatiana Makhalova on 01/07/2020.
//  Copyright © 2020 Tatiana Makhalova. All rights reserved
//
//  The code is a based on the In-Close4 algorithm: Simon Andrews (2016) In-Close4 64 bit version [https://sourceforge.net/projects/inclose/].


#include <iostream>
#include <cstring>  // for memcpy
#include <tuple>
#include <fstream>
#include <iostream>  //for cin >> and cout <<
#include <string>
#include <sstream>
#include <vector>
#include <time.h>  // time_t, struct tm, difftime, time, mktime
#include <mm_malloc.h>
#include <basic_types.h>
#include "tree.hpp"
#include <filesystem>
#include <dirent.h>
#include "io.hpp"
#include <sys/stat.h>

VECTOR_TYPE** A_list;
VECTOR_TYPE** B_list;
int* Iter_list;
int m, n;  // object, attribute number
int mArray, nArray;
std::string fin = "", fout;  // input-output file names
int* rowOriginal;  // maps ham-sorted rows to original order
int* colOriginal;  // maps sorted columns to original order
int* colSup;  // column support (for sorting and skipping empty columns)
std::vector<std::string> attributes;
bool strat_int = true;
bool output_itemsets = true;
int level_id = 0;
int max_level_id = -1;
ALIGNED_VECTOR_TYPE(**context);
ALIGNED_VECTOR_TYPE(**contextTemp);
double* time_splits;
int* level_size;
int* tree_size_array;

void GDPM_error() {  // GDPM requires at least one parameter : full output file name
    exit(0);
}

int main(int argc, char** argv) {
    int  GetNextLevel(int last_cid, int level_id, node* &root_ext);
    int  GetNextLevel_INT(int last_cid, int level_id, node* &root);
    void datFileInput(std::string fin);
    void sortColumnsR();
    void sortColumns();
    void outputConceptsByLevels();
    void sortRows();
    void read_param(int argc, char*argv[], std::string& fin, bool& strat_int,
                    int& max_level_id, bool& output_itemsets);
    std::string strategy_name = "I";  // GDPM-INT
    read_param(argc, argv, fin, strat_int, max_level_id, output_itemsets);
    if (!strat_int) {
        strategy_name = "E";  // GDPM-EXT
    }
    struct stat buffer;
    int check = 0;
    size_t found = fin.find_last_of("//");
    fout  = fin.substr(0, found+1) + "closure_structure/";
    check = 0;
    if (stat(fout.c_str(), &buffer) != 0) {
        check = mkdir(fout.c_str(), 0777);
        std::cout << check << "  " << stat(fout.c_str(), &buffer) << "  " << fout << std::endl;
    }
    if (check == 0) {
        fout += "/dat/";
        check = mkdir(fout.c_str(), 0777);
    }
    if (check != 0) {
        std::cout << "Unable to create an output directory" << std::endl;
    }
    std::cout << "Intput file " << fin << std::endl;
    std::cout << "Output file " << fout << std::endl;
    datFileInput(fin);  // reading dataset
    std::cout << "Dataset of size " << n << " x " << m << std::endl;
    colOriginal = new int[m];
    for (int i = 0; i < m; i++) colOriginal[i] = i;  // init column index array for sorting
    sortColumns();  // sorting columns in asceding order w.r.t. column support
    for (int j = 0; j < m; j++)
        for (int i = 0; i < n; i++)
            if (contextTemp[j][(i>>6)] & (VECTOR_1 << (i % 64)))
                context[i][(j>>6)] |= (VECTOR_1 << (j % 64));
    rowOriginal = new int[n];
    for (int i=0; i < n; i++) rowOriginal[i]=i;
//    sortRows();
    if ((max_level_id == -1) || (max_level_id > m)) {
        max_level_id = m;
    }
    DYN_ALIGNED_DELETE_ARR(contextTemp);
    time_splits = new double[m];
    level_size = new int[m];
    tree_size_array = new int[m];
    A_list = new VECTOR_TYPE* [1];
    B_list = new VECTOR_TYPE* [1];
    node* root = getNewNode();
    VECTOR_TYPE Ac[nArray];
    for (int i = 0; i < nArray; i++) Ac[i] = ~0;
    A_list[0] = &Ac[0];
    VECTOR_TYPE Bc[mArray];
    for (int i = 0; i < mArray; i++) Bc[i] = 0;
    B_list[0] = &Bc[0];
    Iter_list = new int[1];
    Iter_list[0] = 0;
    bool write_zero = false;
    for (int i = 0; i < m; i++) {
        bool add = true;
        int j = 0;
        while ((j < n) && add) {
            if (!(context[j][i>>6] & (VECTOR_1 << (i % 64))))
                add = false;
            j++;
        }
        if (add)
            Bc[i>>6] |= (VECTOR_1 << (i % 64));
        write_zero += add;
    }
    if (write_zero)
        write_itemsets(fout +  std::to_string(0), A_list, B_list, 1, n, m,
                       colOriginal, attributes);
    level_size[0] = 1;
    tree_size_array[0] = 1;
    if (strat_int)
        insert(Bc, root, m);
    else
        insert(Ac, root, n);
    bool run = true;
    int last_cid = 1;
    int total_sum = 0;
    while (run && level_id < max_level_id) {
        int new_last_cid = 0;
        if (strat_int)
            new_last_cid  = GetNextLevel_INT(last_cid, level_id, root);
        else
            new_last_cid  = GetNextLevel(last_cid, level_id, root);
        level_id++;
        total_sum += new_last_cid;
        std::cout << "\n OUTPUT " << level_id << " " << new_last_cid <<
                     "; time  " << time_splits[level_id] << std::endl;
        std::cout << " size : " << n_nodes(root) << std::endl;
        run = new_last_cid > 0;
        if (output_itemsets && run)
                    write_itemsets(fout +  std::to_string(level_id), A_list,
                                   B_list, new_last_cid, n, m, colOriginal,
                                   attributes);
        // read_itemsets("res" +  std::to_string(level_id) + ".csv");
        last_cid = new_last_cid;
        level_size[level_id] = new_last_cid;
        tree_size_array[level_id] = n_nodes(root);
    }
    write_summary(fout + "summary_" + strategy_name + ".csv", time_splits,
                  level_size, tree_size_array, level_id);
}

int  GetNextLevel(int last_cid, int level_id, node* &root_ext) {
    VECTOR_TYPE** A_list_new = new VECTOR_TYPE* [m * (last_cid)];
    VECTOR_TYPE** B_list_new = new VECTOR_TYPE* [m * (last_cid)];
    int* Iter_list_new = new int[m * (last_cid)];
    VECTOR_TYPE* Ac;
    VECTOR_TYPE* Bc;
    int cid_new = 0;
    time_t start, end;
    time(&start);
    for (int cid = 0; cid < last_cid; cid++) {
        Ac = A_list[cid];
        Bc = B_list[cid];
        for (int i = Iter_list[cid]; i < m; i++) {
            if (!(Bc[i>>6] & (VECTOR_1 << (i % 64)))) {  // compute next
                VECTOR_TYPE* Ac_new = new VECTOR_TYPE[nArray];
                for (int j = 0; j < nArray; j++) Ac_new[j] = 0;
                int size = 0;
                for (int j = 0; j < n; j++)
                    if ((Ac[j>>6] & (VECTOR_1 << (j % 64))) &&
                        (context[j][i>>6] & (VECTOR_1 << (i % 64)))) {
                        Ac_new[j>>6] |= (VECTOR_1 << (j % 64));
                        size++;
                    }
                if (size > 0) {
                    node* c_node = insert(Ac_new, root_ext, n);
                    if (c_node != nullptr) {
                        VECTOR_TYPE* Bc_new = new VECTOR_TYPE[mArray];
                        memcpy(Bc_new, Bc, mArray*8);
                        for (int ii = 0; ii < m; ii++) {
                            if (!((Bc[ii>>6] & (VECTOR_1 << (ii % 64))) && (ii != i))) {
                                bool add = true; int jj = 0;
                                while ((jj < n) && add) {
                                    if (Ac_new[jj>>6] & (VECTOR_1 << (jj % 64)))
                                        if (!(context[jj][ii>>6] & (VECTOR_1 << (ii % 64))))
                                            add = false;
                                    jj++;
                                }
                                if (add)
                                    Bc_new[ii>>6] |= (VECTOR_1 << (ii % 64));
                            }
                        }
                        A_list_new[cid_new] = Ac_new;
                        B_list_new[cid_new] = Bc_new;
                        Iter_list_new[cid_new] = i + 1;
                        cid_new++;
                    }
                } else
                    delete [] Ac_new;
           }
        }
    }
    time(&end);
    time_splits[level_id + 1] = double(end - start);
    delete [] A_list;
    delete [] B_list;
    delete [] Iter_list;
    A_list = A_list_new;
    B_list = B_list_new;
    Iter_list = Iter_list_new;
    return cid_new;
}

int  GetNextLevel_INT(int last_cid, int level_id, node* &root) {
    VECTOR_TYPE** A_list_new = new VECTOR_TYPE* [m * (last_cid)];
    VECTOR_TYPE** B_list_new = new VECTOR_TYPE* [m * (last_cid)];
    int* Iter_list_new = new int[m * (last_cid)];
    VECTOR_TYPE* Ac; VECTOR_TYPE* Bc;
    int cid_new = 0;
    time_t start, end;
    time(&start);  /* get current time; same as: timer = time(NULL)  */
    for (int cid = 0; cid < last_cid; cid++) {
        Ac = A_list[cid];
        Bc = B_list[cid];
        for (int i = Iter_list[cid]; i < m; i++) {
            if (!(Bc[i>>6] & (VECTOR_1 << (i % 64)))) {  // compute next
                VECTOR_TYPE* Ac_new = new VECTOR_TYPE[nArray];
                for (int j = 0; j < nArray; j++) Ac_new[j] = 0;
                int size = 0;
                for (int j = 0; j < n; j++)
                    if ((Ac[j>>6] & (VECTOR_1 << (j % 64))) && (context[j][i>>6] & (VECTOR_1 << (i % 64)))) {
                        Ac_new[j>>6] |= (VECTOR_1 << (j % 64));
                        size++;
                    }
                if (size > 0) {
                    VECTOR_TYPE* Bc_new = new VECTOR_TYPE [mArray];
                    memcpy(Bc_new, Bc, mArray*8);
                    for (int ii = 0; ii < m; ii++) {
                        if (!((Bc[ii>>6] & (VECTOR_1 << (ii % 64))) && (ii != i))) {  // all exept of closed and i
                            bool add = true; int jj = 0;
                            while ((jj < n) && add) {
                                if (Ac_new[jj>>6] & (VECTOR_1 << (jj % 64)))
                                    if (!(context[jj][ii>>6] & (VECTOR_1 << (ii % 64))))
                                        add = false;
                                jj++;
                            }
                            if (add)
                                Bc_new[ii>>6] |= (VECTOR_1 << (ii % 64));
                        }
                    }
                    node* c_node = insert(Bc_new, root, m);
                    if (c_node != nullptr) {
                        A_list_new[cid_new] = Ac_new;
                        B_list_new[cid_new] = Bc_new;
                        Iter_list_new[cid_new] = i + 1;
                        cid_new++;
                    } else {
                        delete [] Bc_new;
                    }
                } else {
                    delete [] Ac_new;
                }
           }
        }
    }
    time(&end);
    time_splits[level_id + 1] = double(end - start);
    delete [] A_list;
    delete [] B_list;
    delete [] Iter_list;
    A_list = A_list_new;
    B_list = B_list_new;
    Iter_list = Iter_list_new;
    return cid_new;
}

void datFileInput(std::string fin) {
    std::vector<int> attribute_indices;
    int i, j;  // object and attribute counters
    std::ifstream datFile;
    datFile.open(fin.c_str());
    std::cout << "\nReading data..." << std::endl;;
    n = count(std::istreambuf_iterator<char> (datFile), std::istreambuf_iterator<char> (), '\n');
    nArray = (n - 1)/64 + 1;
    datFile.clear();
    datFile.seekg (0, std::ios::beg);
    std::string line;
    while (getline(datFile, line)) {
        std::istringstream s(line);
        std::string cell;
        while (getline(s, cell, ' ')) {
            std::vector<std::string>::iterator i = find(attributes.begin(), attributes.end(), cell);
            if (i == attributes.end()) {
                attributes.push_back(cell);
                std::istringstream(cell) >> j;
                if (j + 1 > attribute_indices.size()) {
                    while (j > attribute_indices.size()) {
                        attribute_indices.push_back(-1);
                    }
                    attribute_indices.push_back(attributes.size() - 1);
                } else {
                    attribute_indices[j] = attributes.size() - 1;
                }
            }
        }
    }
    m = attributes.size();
    mArray = (m - 1)/64 + 1;
    contextTemp = DYN_ALIGNED_VECTOR_TYPE_PTR_ARR(m);
    for (j = 0; j < m; j++) {
        contextTemp[j] = DYN_ALIGNED_VECTOR_TYPE_ARR(nArray);
        for (i=0; i < nArray; i++) contextTemp[j][i]=0;
    }
    context = DYN_ALIGNED_VECTOR_TYPE_PTR_ARR(n);
    for (i = 0; i < n; i++) {
        context[i] = DYN_ALIGNED_VECTOR_TYPE_ARR(mArray);
        for (j = 0; j < mArray; j++) context[i][j]=0;
    }
    colSup = new int[m];
    for (int i = 0; i < m; i++) colSup[i] = 0;
    datFile.clear();
    datFile.seekg(0, std::ios::beg);

    int r = 0; // row number
    while (getline(datFile, line)) {
        std::istringstream s(line);
        std::string cell;
        while (std::getline(s, cell, ' ')) {
            std::vector<std::string>::iterator k = find(attributes.begin(), attributes.end(), cell);
            if (k != attributes.end()) {
                int j;
                std::istringstream(cell) >> j;
                // set context bit to true where I is byte: i div 8, bit: i mod 8
                contextTemp[attribute_indices[j]][(r>>6)] |= (VECTOR_1 << (r % 64));
                colSup[attribute_indices[j]]++;  // increment column support for attribute j
            }
        }
        r++;
    }
    attribute_indices.erase(attribute_indices.begin(), attribute_indices.end());
    datFile.close();
}
void sortColumns() {  // ascending order
    void colQSort(int colSup[], int colOriginal[], int lo, int high);
    int temp, i, j;
    // bubble sort column indexes (logical sort) - ascending order of support
    for (i = 0 ; i < m ; i++) {
        for (j = 0 ; j <  m-i-1; j++) {
            if (colSup[j] > colSup[j+1]) {
                temp = colSup[j];
                colSup[j] = colSup[j+1];
                colSup[j+1] = temp;
                temp = colOriginal[j];
                colOriginal[j] = colOriginal[j+1];  // keep track of original columns
                colOriginal[j+1] = temp;
            }
        }
    }
    // rewrite sorted context (physical sort)
    int tempColNums[MAX_COLS];
    int rank[MAX_COLS];
    for (j = 0; j < m; j++) {
        tempColNums[j] = colOriginal[j]; // use original col nos to index the sort
        rank[colOriginal[j]] = j;  // record the ranking of the column
    }
    for (j = 0; j < m-1; j++) {
        for (i = 0; i < nArray; i++) {
            VECTOR_TYPE temp = contextTemp[j][i];
            contextTemp[j][i] = contextTemp[tempColNums[j]][i];
            contextTemp[tempColNums[j]][i] = temp;
        }
        tempColNums[rank[j]] = tempColNums[j];  /* make note of where
                                                 swapped-out col has moved to
                                                 using its rank */
        rank[tempColNums[j]] = rank[j];
    }
}

void sortColumnsR() {  // descending order
    void colQSort(int colSup[], int colOriginal[], int lo, int high);
    int temp, i, j;
    // bubble sort column indexes (logical sort) - ascending order of support
    for (i = 0 ; i < m ; i++) {
        for (j = 0; j < m-i-1; j++) {
            if (colSup[j] < colSup[j+1]) {
                temp = colSup[j];
                colSup[j] = colSup[j+1];
                colSup[j+1] = temp;
                temp = colOriginal[j];
                colOriginal[j] = colOriginal[j+1];  // keep track of original columns
                colOriginal[j + 1] = temp;
            }
        }
    }
    // rewrite sorted context (physical sort)
    int* tempColNums = new int[m];
    int* rank = new int[m];
    for (j = 0; j < m; j++) {
        tempColNums[j] = colOriginal[j]; // use original col nos to index the sort
        rank[colOriginal[j]] = j;  // record the ranking of the column
    }
    for (j = 0; j < m-1; j++) {
        for (i = 0; i < nArray; i++) {
            VECTOR_TYPE temp = contextTemp[j][i];
            contextTemp[j][i] = contextTemp[tempColNums[j]][i];
            contextTemp[tempColNums[j]][i] = temp;
        }
        tempColNums[rank[j]] = tempColNums[j];
        rank[tempColNums[j]] = rank[j];
    }
}

void sortRows() {
    void quickBent2Hamsort(int left, int right, int original[]);
    quickBent2Hamsort(0, n-1, rowOriginal);
//    rewrite sorted context (physical sort)
    int* tempRowNums = new int[MAX_ROWS];
    int* rank = new int[MAX_ROWS];
//    for (int i = 0; i < m; i++) {
    for (int i = n - 1; i >= 0; i--) {
        tempRowNums[i] = rowOriginal[i];  // use original row nos to index the sort
        rank[rowOriginal[i]] = i;  // record the ranking of the row
    }
//    for (int i = 0; i < n-1; i++) {
    for (int i = n-1; i >= 0; i--) {
//        for (int j = 0; j < mArray; j++) {
        for (int j = mArray-1; j >= 0; j--) {
            VECTOR_TYPE temp = context[i][j];
            context[i][j] = context[tempRowNums[i]][j];
            context[tempRowNums[i]][j] = temp;
        }
        tempRowNums[rank[i]] = tempRowNums[i];
        rank[tempRowNums[i]] = rank[i];
    }
    delete tempRowNums, rank;
}

void quickBent2Hamsort(int left, int right, int original[]) {
    bool biggerHam(int i1, int i2);
    bool eq(int i1, int i2);
    void exchange(int original[], int from, int to);
    void insertionSort(int left, int right, int original[]);
    int median3(int original[], int i, int j, int k);
    if (right < left) return;
    int n = right - left + 1;
    if (n <= 8) {
        insertionSort(left, right, original);
        return;
    } else if (n <= 60) {
        int m = median3(original, left, left + n/2, right);
        exchange(original, m, left);
    } else {
        int eps = n/8;
        int mid = left + n/2;
        int m1 = median3(original, left, left + eps, left + eps + eps);
        int m2 = median3(original, mid - eps, mid, mid + eps);
        int m3 = median3(original, right - eps - eps, right - eps, right);
        int ninther = median3(original, m1, m2, m3);
        exchange(original, ninther, left);
    }
    int i = left;
    int j = right + 1;
    int p = left;
    int q = right + 1;
    int pivotIndex = left;
    int pivotValue =  original[pivotIndex];
    while (true) {
        while (biggerHam(original[++i], pivotValue))
            if (i == right)
                break;
        while (biggerHam(pivotValue, original[--j]))
            if (j == left)
                break;
        // pointers cross
        if (i == j && eq(original[i], pivotValue))
            exchange(original, ++p, i);
        if (i >= j)
            break;
        exchange(original, i, j);
        if (eq(original[i], pivotValue))
            exchange(original, ++p, i);
        if (eq(original[j], pivotValue))
            exchange(original, --q, j);
    }
    // exhange equal parts on left and right sides to the middle
    i = j + 1;
    for (int k = left; k <= p; k++)
        exchange(original, k, j--);
    for (int k = right; k >= q; k--)
        exchange(original, i++, k);
    quickBent2Hamsort(left, j, original);
    quickBent2Hamsort(i, right, original);
    }

void insertionSort(int left, int right, int original[]) {
    bool biggerHam(int i1, int i2);
    void exchange(int original[], int left, int right);
    for (int i = left; i <= right; i++) {
        for (int j = i;  j > left && biggerHam(original[j], original[j-1]); j--)
            exchange(original, j, j - 1);
    }
}

int median3(int original[], int i, int j, int k) {
    bool biggerHam(int i1, int i2);
    int med;
    biggerHam(original[i], original[j]) ?
               (biggerHam(original[j], original[k]) ? med = j :
                biggerHam(original[i], original[k]) ? med = k : med = i) :
               (biggerHam(original[k], original[j]) ? med = j :
                biggerHam(original[k], original[i]) ? med = k : med = i);
    return med;
}

bool biggerHam(int i1, int i2) {
    // for (int j=0;j<nArray;j++) {
    for (int j=nArray-1; j >= 0; j--) {
        if (context[i1][j] > context[i2][j]) return true;
        if (context[i2][j] > context[i1][j]) return false;
    }
    return false;
}

bool eq(int i1, int i2) {
    for (int j=nArray-1; j >= 0; j--) {
        if (context[i1][j] != context[i2][j])
            return false;
    }
    return true;
}

void exchange(int original[], int from, int to) {
    int temp = original[from];
    original[from] = original[to];
    original[to] = temp;
}

void read_param(int argc, char* argv[], std::string& fin, bool& strat_int,
                int& max_level_id, bool& output_itemsets) {
    if ( argc < 2 ) { GDPM_error(); return; }
    fin = argv[1];
    std::cout << argc << "\n";
    int c = 2;
    if (c < argc) {
        if ( strchr (argv[c], 'E') ) strat_int = false;
        if ( strchr (argv[c], 'Q') ) output_itemsets = false;
        c++;
        while (c < argc) {
            switch (argv[c][0]) {
                case 'E': strat_int = false; break;
                case 'Q': output_itemsets = false; break;
                case '-':  if (argv[c][1] == 'n') {
                    max_level_id = std::stoi(argv[c + 1]);
                    c++;
                } break;
            }
            c++;
        }
    }
}
