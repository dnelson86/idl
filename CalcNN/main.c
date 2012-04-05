//
// CalcNN (nearest neighbor) Routine for IDL
// dnelson jan.2012
// based on CalcHsml (for Python) by Mark Vogelsberger
//

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "idl_export.h"
#include "main.h"

// main routine in C
int CalcNN_natural(int NumPart, int NumOrig, float* data_pos, float* src_pos, int* ind_out)
{
  int i, ind, nNodes;
  float pos[3];
  MyDouble mindist;
  int guess = -1;

  P = (struct particle_data *) malloc (NumPart * sizeof (struct particle_data));
	
  // copy data into P
  for (i=0; i < NumPart; i++)
  {
    P[i].Pos[0] = data_pos[3*i+0];
    P[i].Pos[1] = data_pos[3*i+1];
    P[i].Pos[2] = data_pos[3*i+2];
    P[i].Mass   = 1.0;
    P[i].Type   = 0; // search for type 0 (gas)
  }

  // allocate and build tree
  tree_treeallocate(4.0 * NumPart, NumPart);
  nNodes = tree_treebuild();
	
#ifdef VERBOSE
  printf("Built tree, nodes = %d\n", nNodes);

  // use tree for neighbor searches
  printf("finding neighbours...\n [");
  
  int signal = 0;
#endif
	
  for (i=0; i < NumOrig; i++)
  {
#ifdef VERBOSE
    // output progress marker
    if (i > (signal / 100.0) * NumOrig)
    {
      printf("x");
      fflush(stdout);
      signal++;
    }
#endif
    // set search position and find nearest neighbor
    pos[0] = src_pos[3*i+0];
    pos[1] = src_pos[3*i+1];
    pos[2] = src_pos[3*i+2];
    mindist = -1;

    ind = ngb_treefind_nearest_local(pos,-1,&mindist,guess);		

    // store result
    ind_out[i] = ind;

    // due to peano hilbert ordering, set last result as guess for next (much faster)
    guess = ind;
  }

  tree_treefree();
  free(P);

  IF_VERBOSE(printf("]\ndone.\n"));
	
  return 1;
}

/* IDL glue routine
 
  example:
    
    ; prepare inputs
    NumPart = long(nSrcTargs)
    nSrcOrigs = long(nSrcOrigs)

    Pos     = fltarr(3,nSrcTargs)
    PosSrc  = fltarr(3,nSrcOrigs)

    BoxSize = float(10000.0) ;kpc
    
    ind_out = lonarr(NumPart)
    
    ; call CalcNN
    ret = Call_External('/n/home07/dnelson/idl/CalcNN/CalcNN.so', 'CalcNN', $
                        NumPart,nSrcOrigs,Pos,PosSrc,BoxSize,ind_out,/CDECL)
*/

int CalcNN(int argc, void* argv[])
{
  float *Pos,*PosSrc;
  int *ind_out;
  int NumSrc;
	
  char buf[128];

  // validate input
  if (argc != 6)
  {
    sprintf(buf,"Wrong number of arguments (%d)!\n",argc);
    IDL_Message(IDL_M_GENERIC,IDL_MSG_RET,buf);
    return 0;
  } else {
#ifdef TWODIMS
    IDL_Message(IDL_M_GENERIC,IDL_MSG_RET,"CalcNN Loaded (TWODIMS!).");
#else
    IDL_Message(IDL_M_GENERIC,IDL_MSG_RET,"CalcNN Loaded (NDIMS=3!).");
#endif
  }

  // should check: pos has 2 dims, first dim is 3
  // should check: mass has 1 dim
	
  // inputs and return by reference
  NumPart = *(int *)argv[0];
  NumSrc  = *(int *)argv[1];
  Pos     = (float *)argv[2];
  PosSrc  = (float *)argv[3];
  BoxSize = *(float *)argv[4];
  ind_out = (int *)argv[5];

  // defaults
  Softening = 1.0;
  DesNumNgb = 32;
  DesNumNgbDev = 0;

#ifdef VERBOSE
  printf("Input data:\n");
  printf("NumPart      = %d\n", NumPart);
  printf("NumSrc       = %d\n", NumSrc);
  printf("BoxSize      = %g\n", BoxSize);
#endif

  // calculate NNs
  return CalcNN_natural(NumPart,NumSrc,Pos,PosSrc,ind_out);
}
