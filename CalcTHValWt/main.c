/*
 * CalcTHValWt Routine for IDL
 * dnelson may.2012
 * same as CalcTHVal but accepts weights for each value (cell masses for Arepo)
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "idl_export.h"
#include "main.h"

// main routine in C
int CalcTHValWt_natural(int NumPart, float* data_pos_val_wt, int NumSearch, float* search_pos, 
                        int DesNumNgb, int DesNumngbDev, float BoxSize, float* val_out)
{
  int i;
  int nNodes;
  float hsml = 1.0, val;
  float pos[3];

  // instead of copying, just cast the input data as an array of particle_data structs
  P = (struct particle_data *) data_pos_val_wt;

  // allocate and build tree
  tree_treeallocate(2.0 * NumPart, NumPart);
  nNodes = tree_treebuild();

#ifdef VERBOSE
  // use tree for neighbor searches
  printf("Built tree with %d nodes. Finding neighbors...\n [", nNodes);

  int signal = 0;
#endif
	
  for (i=0; i < NumSearch; i++)
  {
    val = 0.0;
#ifdef VERBOSE
    // output progress marker
    if (i > (signal / 100.0) * NumSearch)
    {
      printf("x");
      fflush(stdout);
      signal++;
    }
#endif
    // get current search position
    pos[0] = search_pos[3*i+0];
    pos[1] = search_pos[3*i+1];
    pos[2] = search_pos[3*i+2];

    // calculate smoothing length and use this value as a guess for the next point
    hsml = ngb_treefind(pos, hsml * 1.0, &val);

    // save value for this search position
    val_out[i] = val;
  }

  tree_treefree();

  IF_VERBOSE(printf("]\n"));
	
  return 1;
}

/* IDL glue routine
 
  example:

    nNGB = 32
    
    ; prepare inputs
    NumPart   = long(npts)
    PosVal    = fltarr(5,npts) ; [x,y,z,function value,element weight (mass)]
    NumSrc    = long(nsearchpts)
    SearchPos = fltarr(3,nsearchpts)
    DesNumNgb    = long(nNGB)
    BoxSize      = 0.0
    TophatMode   = 1
    
    val_out = fltarr(NumSrc)
    
    ; call CalcTHVal
    ret = Call_External('/n/home07/dnelson/idl/CalcTHValWt/CalcTHValWt_3D.so', 'CalcTHVal', $
                        NumPart,PosVal,NumSrc,SearchPos,DesNumNgb,TophatMode,BoxSize,val_out,/CDECL)
*/

int CalcTHValWt(int argc, void* argv[])
{
  int NumSearch;
  float *PosValWt,*SrcPos,*val_out;
	
  char buf[128];

  // validate input
  if (argc != 8)
  {
    sprintf(buf,"Wrong number of arguments (%d)!\n",argc);
    IDL_Message(IDL_M_GENERIC,IDL_MSG_RET,buf);
    return 0;
  } else {
#ifdef   ONEDIMS
IDL_Message(IDL_M_GENERIC,IDL_MSG_RET,"CalcTHValWt Loaded (ONEDIMS!).");
#else
#ifndef  TWODIMS
IDL_Message(IDL_M_GENERIC,IDL_MSG_RET,"CalcTHValWt Loaded (NDIMS=3!).");
#else
IDL_Message(IDL_M_GENERIC,IDL_MSG_RET,"CalcTHValWt Loaded (TWODIMS!).");
#endif
#endif /* ONEDIMS */

  }

  // should check: pos has 2 dims, first dim is 3	
  // inputs and return by reference
  NumPart   = *(int *)argv[0];
  PosValWt  = (float *)argv[1];
  NumSearch = *(int *)argv[2];
  SrcPos    = (float *)argv[3];
	
  DesNumNgb    = *(int *)argv[4];
  DesNumNgbDev = 0; // no deviation allowed
  TophatMode   = *(int *)argv[5]; // 1=mean, 2=total, 3=total/volume (density)
  BoxSize      = *(float *)argv[6];
  Softening    = 1.0; // code units

  val_out = (float *)argv[7];

  // allocate Ngblist (up to NumPart big?)
  Ngblist = (int *) malloc(NumPart * sizeof(int));

#ifdef VERBOSE
  //printf("Input data:\n");
  //printf("NumPart      = %d\n", NumPart);
  //printf("NumSearch    = %d\n", NumSearch);
  //printf("DesNumNgb    = %d\n", DesNumNgb);
  //printf("DesNumNgbDev = %d\n", DesNumNgbDev);
  //printf("BoxSize      = %g\n", BoxSize);
#endif

  // calculate tophat values
  int status = CalcTHValWt_natural(NumPart,PosValWt,NumSearch,SrcPos,DesNumNgb,DesNumNgbDev,BoxSize,val_out);

  // free Ngblist and return
  //free(Ngblist);
  return status;
}
