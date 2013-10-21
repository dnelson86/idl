#ifndef ALLVARS_H
#define ALLVARS_H

#define VERBOSE

#ifdef VERBOSE
#define IF_VERBOSE(expr) (expr)
#else
#define IF_VERBOSE(expr) ((void)0)
#endif

// single or double precision
#ifdef INT64_PRECISION
typedef long long MyInt;
typedef int MyIndType; // can change if we need to sort more than max(int) num elements
#else
#ifdef UINT64_PRECISION
typedef unsigned long long MyInt;
typedef int MyIndType;
#else
// INT32_PRECISION
typedef int MyInt;
typedef int MyIndType;
#endif
#endif

// globals
MyIndType NumData;
MyInt *DataIn;
MyIndType *inds_out;
MyInt *A, *B;
	
// data,index pairs
struct sort_data
{
  MyInt ID;
  MyIndType index;	
} *data;

// definitions
typedef int (*comparer)(const void *, const void *);
void* my_qsort (void *arr, size_t length, size_t element_size, comparer cmp);
void _qsort_SD (void *arr, int left, int right, comparer cmp);

#endif
