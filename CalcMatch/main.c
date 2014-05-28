/*
 * CalcMatch Routine for IDL
 * dnelson feb.2014
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
  //return ( id1 - id2 ); // if id1-id2>(32 bits) the overflow kills the return!
  if( id1 > id2 ) return 1;
  if( id1 < id2 ) return -1;
  return 0;
}

/*
int inplace_indcomp_func(const void *a, const void *b) // unused
{
  MyIndType ind1 = *(MyIndType *)a;
  MyIndType ind2 = *(MyIndType *)b;
  return ( ind1 - ind2 );
}
*/

// sort the index array, and in the comparator use it to index the actual IDs

int indexing_comp_func(const void *a, const void *b)
{
  MyInt id1 = DataIn[ *(MyIndType *)a ];
  MyInt id2 = DataIn[ *(MyIndType *)b ];
  //return ( id1 - id2 );
  if( id1 > id2 ) return 1;
  if( id1 < id2 ) return -1;
  return 0;
}

int indexA_comp_func(const void *a, const void *b)
{
  MyInt id1 = A[ *(MyIndType *)a ];
  MyInt id2 = A[ *(MyIndType *)b ];
  //return ( id1 - id2 );
  if( id1 > id2 ) return 1;
  if( id1 < id2 ) return -1;
  return 0;
}
int indexB_comp_func(const void *a, const void *b)
{
  MyInt id1 = B[ *(MyIndType *)a ];
  MyInt id2 = B[ *(MyIndType *)b ];
  //return ( id1 - id2 );
  if( id1 > id2 ) return 1;
  if( id1 < id2 ) return -1;
  return 0;
}

// variant of indexB which uses the addresses of B to insure stability for the CalcMatchDupe case

int indexB_comp_stable_func(const void *a, const void *b)
{
  MyInt id1 = B[ *(MyIndType *)a ];
  MyInt id2 = B[ *(MyIndType *)b ];
  
  if ( id1 != id2 ) {
    if( id1 > id2 ) return 1;
    return -1;
    //return ( id1 - id2 ); // data value difference
  }
  
  //return &( B[ *(MyIndType *)a ] ) - &( B[ *(MyIndType *)b ] ); // address difference
  
  if( &(B[ *(MyIndType *)a ]) > &(B[ *(MyIndType *)b ]) ) return 1;
  return -1;
}

/* 64bit ID with 32bit indices sort example, return permutation indices:

    ; prepare inputs
    NumData  = long64(npts)
    data     = lon64arr(npts) ; particle/tracer IDs
    inds_out = lindgen(npts) ; array indices, input as sequential
    method   = 3L ; 3=fhtr sort_data specialized
    
    ret = Call_External('/n/home07/dnelson/idl/CalcMatch/CalcMatch_int64.so', 'CalcSort', $
                        NumData,data,inds_out,method,/CDECL)

  64bit ID with 32bit indices sort example, modify inplace

    ; prepare inputs
    NumData  = long64(npts)
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
  
  // inputs and return by reference
  NumData  = *(long long *)argv[0];

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
    numA  = long64(n_elements(A))
    numB  = long64(n_elements(B))
    A     = long64(A) ; particle/tracer IDs
    B     = long64(B) ; particle/tracer IDs
    
    inds_A_out = lon64arr(n_elements(A))
    inds_B_out = lon64arr(n_elements(B))
    
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

/* match two 64bit ID arrays, where B can have duplicates (i.e. B is tracer_parentIDs, A is gasIDs)

    ; prepare inputs/outputs
    numA  = long64(n_elements(A))
    numB  = long64(n_elements(B))
    
    inds_A_out = l64indgen(numA)
    inds_B_out = l64indgen(numB)
    count      = -1L
    
    ret = Call_External('/n/home07/dnelson/idl/CalcMatch/CalcMatch_int64.so', 'CalcMatchDupe', $
                        numA,numB,A,B,inds_A_out,inds_B_out,count,/CDECL)
    
    ; take index subsets
    inds_A_out = inds_A_out            ; child_counts for each parent (number of B duplicates per A)
    inds_B_out = inds_B_out[0:count-1] ; indices of B matched to A, ordered in their original orders

*/

int CalcMatchDupe(int argc, void* argv[])
{
  MyIndType i = 0, j = 0;
  char buf[128];
  
  MyIndType numA, numB;
  MyIndType *inds_A_out, *inds_B_out, *count_out;
  MyIndType *inds_A, *inds_B;

  // validate input
  if (argc != 7)
  {
    sprintf(buf,"Wrong number of arguments (%d)!\n",argc);
    IDL_Message(IDL_M_GENERIC,IDL_MSG_RET,buf);
    return 0;
  }
  
  // inputs and return by reference
  numA   = *(MyIndType *)argv[0];
  numB   = *(MyIndType *)argv[1];
  
  A = (MyInt *)argv[2];
  B = (MyInt *)argv[3];
  
  inds_A_out = (MyIndType *)argv[4];
  inds_B_out = (MyIndType *)argv[5];
  count_out  = (MyIndType *)argv[6];
   
  // make a copy of the input indices (which are 0,1,2,...)
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
  my_qsort( inds_A, numA, sizeof(MyIndType), indexA_comp_func );

  my_qsort( inds_B, numB, sizeof(MyIndType), indexB_comp_stable_func );

  // allocate space for the offset table
  MyIndType *offsetTable = (MyIndType *) malloc( (numA+1) * sizeof(MyIndType) );
  offsetTable[0] = 0;

  // zero the input/output inds_A_out array so child counts start at zero
  memset( inds_A_out, 0, numA * sizeof(MyIndType) );

  MyIndType indOffset = 0; // sequential walk through inds_B_out
  MyIndType count = 0; // number of children for this parent

  // single walk match, B can contain duplicates
  while (i < numA && j < numB)
  {
    if ( A[ inds_A[i] ] == B[ inds_B[j] ] ) {
      // matched element
      inds_A_out[inds_A[i]] += 1; // child_count of this A parent (in -unsorted- order)
      inds_B_out[indOffset] = inds_B[j]; // ordered child tracer indices
      
      // move forward only in B, in case we have another match with this same A element
      j += 1;
      indOffset += 1;
      count += 1;
    } else if ( A[ inds_A[i] ] < B[ inds_B[j] ] ) {
      // element mismatch, forward in A, reset local child count, update offsetTable
      i += 1;
      offsetTable[i] = offsetTable[i-1] + count;
      count = 0;
      //if ( i % 100 == 0 ) { printf("i=%d j=%d\n",i,j); }
    } else {
      // element mismatch, forward in B
      j += 1;
    }
  }
  
  // number of total matched children
  count_out[0] = indOffset;
  
  // reverse the A sort
  MyIndType *inds_A_back = (MyIndType *) malloc( numA * sizeof(MyIndType) );
    
  for( i = 0; i < numA; i++)
    inds_A_back[inds_A[i]] = i;
  
  // need to rearrange inds_B_out based on the sorting of A (using the offsetTable)
  // first copy the first [0:indOffset-1] of inds_B_out over inds_B (no longer needed)
  memcpy( inds_B, inds_B_out, indOffset * sizeof(MyIndType) );
  
  // walk through A (in original order) and copy its matched B indices sequentially into inds_B_out
  indOffset = 0;
  
  MyIndType curInd;
  int numMatches;
  
  for ( i = 0; i < numA; i++ )
  {
    curInd = inds_A_back[i];
    numMatches = inds_A_out[i];
    
    if ( !numMatches )
      continue;
      
    for ( j = 0; j < numMatches; j++ )
      inds_B_out[indOffset + j] = inds_B[offsetTable[curInd] + j];

    indOffset += numMatches;
  }
  
  free(inds_A_back);
  free(inds_A);
  free(inds_B);

  return 1;
}
