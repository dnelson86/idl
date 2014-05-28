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
typedef long long MyIndType;
#else
#ifdef UINT64_PRECISION
typedef unsigned long long MyInt;
typedef long long MyIndType;
#else
#ifdef UINT32_PRECISION
typedef unsigned int MyInt;
typedef int MyIndType;
#else
// INT32_PRECISION
typedef int MyInt;
typedef int MyIndType;
#endif
#endif
#endif

// globals
long long NumData;
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
void* my_qsort (void *arr, MyIndType length, size_t element_size, comparer cmp);
void _qsort_SD (void *arr, MyIndType left, MyIndType right, comparer cmp);

#endif
