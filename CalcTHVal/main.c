/*
 * CalcTHVal Routine for IDL
 * dnelson apr.2012
 * similar to CalcHSMLds in that it uses a tophat filter of constant neighbor number on a "different set"
 * of points in R^3 than the particles themselves. instead of returning the hsml of the tophat this routine
 * estimates the (mean) value of some quantity specified for each particle, e.g. temperature
 * 1D,2D,3D versions, for periodic set BoxSize>0 (LONG_XYZ not supported though)
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "idl_export.h"
#include "main.h"

// main routine in C
int CalcTHVal_natural(int NumPart, float* data_pos_val, int NumSearch, float* search_pos, 
                      int DesNumNgb, int DesNumngbDev, float BoxSize, float* val_out)
{
  int i;
  int nNodes;
  float hsml = 1.0, val = 0.0;
  float pos[3];

  // instead of copying, just cast the input data as an array of particle_data structs
  P = (struct particle_data *) data_pos_val;

  //for (i=0; i < NumPart; i++)
  //  printf("[%d] %g %g %g val %g\n",i,P[i].Pos[0],P[i].Pos[1],P[i].Pos[2],P[i].Value);

  // allocate and build tree
  tree_treeallocate(4.0 * NumPart, NumPart);
  nNodes = tree_treebuild();
	
#ifdef VERBOSE
  // use tree for neighbor searches
  printf("Built tree with %d nodes. Finding neighbors...\n [", nNodes);

  int signal = 0;
#endif
	
  for (i=0; i < NumSearch; i++)
  {
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
    hsml = ngb_treefind(pos, hsml * 1.1, &val);

    // save value for this search position
    val_out[i] = val;
    //printf("search [%d] %g %g %g hsml %g val %g\n",i,pos[0],pos[1],pos[2],hsml,val);
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
    PosVal    = fltarr(4,npts) ; [x,y,z,function value]
    NumSrc    = long(nsearchpts)
    SearchPos = fltarr(3,nsearchpts)
    DesNumNgb    = long(nNGB)
    DesNumNgbDev = long(0)
    BoxSize      = 0.0
    
    val_out = fltarr(NumSrc)
    
    ; call CalcTHVal
    ret = Call_External('/n/home07/dnelson/idl/CalcTHVal/CalcTHVal_3D.so', 'CalcTHVal', $
                        NumPart,PosVal,NumSrc,SearchPos,DesNumNgb,DesNumNgbDev,BoxSize,val_out,/CDECL)
*/

int CalcTHVal(int argc, void* argv[])
{
  int NumSearch;
  float *PosVal,*SrcPos,*val_out;
	
  char buf[128];

  // validate input
  if (argc != 8)
  {
    sprintf(buf,"Wrong number of arguments (%d)!\n",argc);
    IDL_Message(IDL_M_GENERIC,IDL_MSG_RET,buf);
    return 0;
  } else {
#ifdef   ONEDIMS
IDL_Message(IDL_M_GENERIC,IDL_MSG_RET,"CalcHSMLds Loaded (ONEDIMS!).");
#else
#ifndef  TWODIMS
IDL_Message(IDL_M_GENERIC,IDL_MSG_RET,"CalcHSMLds Loaded (NDIMS=3!).");
#else
IDL_Message(IDL_M_GENERIC,IDL_MSG_RET,"CalcHSMLds Loaded (TWODIMS!).");
#endif
#endif /* ONEDIMS */

  }

  // should check: pos has 2 dims, first dim is 3	
  // inputs and return by reference
  NumPart   = *(int *)argv[0];
  PosVal    = (float *)argv[1];
  NumSearch = *(int *)argv[2];
  SrcPos    = (float *)argv[3];
	
  DesNumNgb    = *(int *)argv[4];
  DesNumNgbDev = *(int *)argv[5];
  BoxSize      = *(float *)argv[6];
  Softening    = 1.0; // code units
	
  val_out = (float *)argv[7];

  // allocate Ngblist
  Ngblist = (int *) malloc((DesNumNgb+DesNumNgbDev) * sizeof(int));

#ifdef VERBOSE
  //printf("Input data:\n");
  //printf("NumPart      = %d\n", NumPart);
  //printf("NumSearch    = %d\n", NumSearch);
  //printf("DesNumNgb    = %d\n", DesNumNgb);
  //printf("DesNumNgbDev = %d\n", DesNumNgbDev);
  //printf("BoxSize      = %g\n", BoxSize);
#endif

  // calculate tophat values
  int status = CalcTHVal_natural(NumPart,PosVal,NumSearch,SrcPos,DesNumNgb,DesNumNgbDev,BoxSize,val_out);

  // free Ngblist and return
  free(Ngblist);
  return status;
}