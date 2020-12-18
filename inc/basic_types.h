

#ifndef INC_BASIC_TYPES_H_
#define INC_BASIC_TYPES_H_
#include <mm_malloc.h>
#define VECTOR_TYPE unsigned long long
#define VECTOR_1  1LL
#define VECTOR_0  ~1LL
#define ALIGNMENT __attribute__ ((aligned (32)))
#define ALIGNED_VECTOR_TYPE(n) VECTOR_TYPE n ALIGNMENT
#define MAX_ATTRIBUTES 40000000        // max number of attributes

#define DYN_ALIGNED_VECTOR_TYPE_ARR(n) reinterpret_cast<VECTOR_TYPE *> (_mm_malloc(n*sizeof(VECTOR_TYPE), 32))
#define DYN_ALIGNED_VECTOR_TYPE_PTR_ARR(n) reinterpret_cast<VECTOR_TYPE **> (_mm_malloc(n*sizeof(VECTOR_TYPE *), 32))
#define DYN_ALIGNED_DELETE_ARR(n) _mm_free(n)


#define MAX_CONS 23000000       // max number of concepts
#define MAX_COLS 5000           // max number of attributes
#define MAX_ROWS 200000         // max number of objects
#define MAX_FOR_B 40000000      // memory for storing intents
#define MAX_FOR_A 300000000     // memory for storing extents

#endif   // INC_BASIC_TYPES_H_
