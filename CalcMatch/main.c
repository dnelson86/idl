/*
 * CalcMatch Routine for IDL
 * dnelson mar.2013
 * replicate behavior of IDL sort() and match() but hopefully faster
 * specialized to 32 and 64 bit integer types (IDs)
 * NOTE: currently uses int for indexing, so restricted to arrays of size <~ 2 bil
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "idl_export.h"
#include "main.h"

// comparison functions

int inplace_comparison_func(const void *a, const void *b)
{
  MyInt id1 = *(MyInt *)a;
  MyInt id2 = *(MyInt *)b;
  return ( id1 - id2 );
}

// sort the index array, and in the comparator use it to index the actual IDs

int indexing_comp_func(const void *a, const void *b)
{
  MyInt id1 = DataIn[ *(MyIndType *)a ];
  MyInt id2 = DataIn[ *(MyIndType *)b ];
  return ( id1 - id2 );
}

int indexA_comp_func(const void *a, const void *b)
{
  MyInt id1 = A[ *(MyIndType *)a ];
  MyInt id2 = A[ *(MyIndType *)b ];
  return ( id1 - id2 );
}
int indexB_comp_func(const void *a, const void *b)
{
  MyInt id1 = B[ *(MyIndType *)a ];
  MyInt id2 = B[ *(MyIndType *)b ];
  return ( id1 - id2 );
}

/* 64bit ID with 32bit indices sort example, return permutation indices:

    ; prepare inputs
    NumData  = long(npts)
    data     = lon64arr(npts) ; particle/tracer IDs
    inds_out = lindgen(npts) ; array indices, input as sequential
    method   = 3L ; 3=fhtr sort_data specialized
    
    ret = Call_External('/n/home07/dnelson/idl/CalcMatch/CalcMatch_int64.so', 'CalcSort', $
                        NumData,data,inds_out,method,/CDECL)

  64bit ID with 32bit indices sort example, modify inplace

    ; prepare inputs
    NumData  = long(npts)
    data     = lon64arr(npts)
    inds_out = [0]
    method   = 12L ; 10=inplace + 2=fhtr (12)

    ret = Call_External('/n/home07/dnelson/idl/CalcMatch/CalcMatch_int64.so', 'CalcSort', $
                        NumData,data,inds_out,method,/CDECL)

*/

int CalcSort(int argc, void* argv[])
{
  int method = 0;
  int inPlace = 0;
  char buf[128];

  // validate input
  if (argc != 4)
  {
    sprintf(buf,"Wrong number of arguments (%d)!\n",argc);
    IDL_Message(IDL_M_GENERIC,IDL_MSG_RET,buf);
    return 0;
  }
  
#ifdef VERBOSE
#ifdef INT64_PRECISION
  //IDL_Message(IDL_M_GENERIC,IDL_MSG_RET,"CalcSort Loaded (INT64).");
#else
  //IDL_Message(IDL_M_GENERIC,IDL_MSG_RET,"CalcSort Loaded (INT32).");
#endif
#endif

  // inputs and return by reference
  NumData  = *(MyIndType *)argv[0];

  // (double memory usage)
  DataIn   = (MyInt *)argv[1];
  inds_out = (MyIndType *)argv[2];
  method   = *(int *)argv[3];  

  if(method <= 0) {
    sprintf(buf,"Error: Method zero or unspecified.");
    IDL_Message(IDL_M_GENERIC,IDL_MSG_RET,buf);
    return 0;
  }

  // if method>10 that means do the sort inplace, return sorted array and not indices
  if(method > 10) {
    inPlace = 1;
    method -= 10;
  }

  // INDICES, INPUT UNMODIFIED, INDS_OUT TO CONTAIN PERMUTATION
  if (!inPlace)
  {
    // glibc sort
    if(method == 1)
      qsort( inds_out, NumData, sizeof(MyIndType), indexing_comp_func );
      //qsort( data, NumData, sizeof(struct sort_data), indexing_comp_func );

    // fhtr sort (threaded, generic)
    if(method == 2)
      my_qsort( inds_out, NumData, sizeof(MyIndType), indexing_comp_func );
  }
  else
  { // INPLACE, REPLACE INPUT WITH SORTED ARRAY
    // glibc sort
    if(method == 1)
      qsort( DataIn, NumData, sizeof(MyInt), inplace_comparison_func );

    // fhtr sort (threaded, generic or specialized to 4/8 byte MyInt)
    if(method == 2 || method == 3)
      my_qsort( DataIn, NumData, sizeof(MyInt), inplace_comparison_func );
  }

  return 1;
}

/* match two 64bit ID arrays, return permutation indices, identical behavior to match.pro:

    ; prepare inputs
    method = 0L ; unused
    numA  = long(n_elements(A))
    numB  = long(n_elements(B))
    A     = long64(A) ; particle/tracer IDs
    B     = long64(B) ; particle/tracer IDs
    
    inds_A_out = lonarr(n_elements(A))
    inds_B_out = lonarr(n_elements(B))
    
    ret = Call_External('/n/home07/dnelson/idl/CalcMatch/CalcMatch_int64.so', 'CalcMatch', $
                        method,numA,numB,A,B,inds_A_out,inds_B_out,/CDECL)
    
    ; take index subsets (note: numA=numB=numMatch replaced in CalcMatch)
    inds_A_out = inds_A_out[0:numA-1]
    inds_B_out = inds_B_out[0:numB-1]

*/

int CalcMatch(int argc, void* argv[])
{
  MyIndType i = 0, j = 0;
  MyIndType indOffset = 0;
  int method = 0;
  char buf[128];
  
  MyIndType numA, numB;
  MyIndType *inds_A_out, *inds_B_out, *count_out;
  MyIndType *inds_A, *inds_B;

  // validate input
  if (argc != 8)
  {
    sprintf(buf,"Wrong number of arguments (%d)!\n",argc);
    IDL_Message(IDL_M_GENERIC,IDL_MSG_RET,buf);
    return 0;
  }
  
#ifdef VERBOSE
#ifdef INT64_PRECISION
  //IDL_Message(IDL_M_GENERIC,IDL_MSG_RET,"CalcSort Loaded (INT64).");
#else
  //IDL_Message(IDL_M_GENERIC,IDL_MSG_RET,"CalcSort Loaded (INT32).");
#endif
#endif

  // inputs and return by reference
  method = *(int *)argv[0]; 
  numA   = *(MyIndType *)argv[1];
  numB   = *(MyIndType *)argv[2];
  
  A = (MyInt *)argv[3];
  B = (MyInt *)argv[4];
  
  inds_A_out = (MyIndType *)argv[5];
  inds_B_out = (MyIndType *)argv[6];
  
  count_out = (MyIndType *)argv[7];
   
  if(method <= 0) {
    sprintf(buf,"Error: Method zero or unspecified.");
    IDL_Message(IDL_M_GENERIC,IDL_MSG_RET,buf);
    return 0;
  }
  
  // make a copy of the indices and sort those
  inds_A = (MyIndType *) malloc(numA * sizeof(MyIndType));
  inds_B = (MyIndType *) malloc(numB * sizeof(MyIndType));
  
  if ( inds_A == NULL || inds_B == NULL ) {
    sprintf(buf,"Error: Failed to allocate memory for indices copy.");
    IDL_Message(IDL_M_GENERIC,IDL_MSG_RET,buf);
    return 0;
  }
  
  memcpy( inds_A, inds_A_out, numA * sizeof(MyIndType) );
  memcpy( inds_B, inds_B_out, numB * sizeof(MyIndType) );
  
  // fhtr sort (threaded), skip for A if method >= 10
  if(method < 10)
    my_qsort( inds_A, numA, sizeof(MyIndType), indexA_comp_func );
    
  my_qsort( inds_B, numB, sizeof(MyIndType), indexB_comp_func );

  // single walk match
  while (i < numA && j < numB)
  {
    if ( A[ inds_A[i] ] == B[ inds_B[j] ] ) {
      // matched element
      inds_A_out[indOffset] = inds_A[i];
      inds_B_out[indOffset] = inds_B[j];
      i += 1; j += 1;
      indOffset += 1;
    } else if ( A[ inds_A[i] ] < B[ inds_B[j] ] ) {
      i += 1;
    } else {
      j += 1;
    }
  }
  
  free(inds_A);
  free(inds_B);
  
  // copy number of matched elements over numA and numB (note: off1=off2)
  count_out[0] = indOffset;

  return 1;
}
