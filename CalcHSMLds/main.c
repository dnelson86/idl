/*
 * CalcHSMLds Routine for IDL
 * dnelson apr.2012
 * unlike CalcHSML this esimates the smoothing lengths (e.g. tophat filter sizes) for a "different set"
 * of points in R^3 which are in general different from the particle positions
 * 1D,2D,3D versions, for periodic set BoxSize>0 (LONG_XYZ not supported though)
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "idl_export.h"
#include "main.h"

// main routine in C
int CalcHSMLds_natural(int NumPart, float* data_pos, int NumSearch, float* search_pos, 
                     int DesNumNgb, int DesNumngbDev, float BoxSize, float* hsml_out)
{
  int i;
  int nNodes;
  float hsml = 1.0;
  float pos[3];

/*
 * P = (struct particle_data *) malloc (NumPart * sizeof (struct particle_data));	
 * // copy data into P
 * for (i=0; i < NumPart; i++)
 * {
 *   P[i].Pos[0] = data_pos[3*i+0];
 *   P[i].Pos[1] = data_pos[3*i+1];
 *   P[i].Pos[2] = data_pos[3*i+2];
 *   //P[i].Mass   = data_mass[i];
 * }
 */

  // instead of copying, just cast the input data as an array of particle_data structs
  P = (struct particle_data *) data_pos;

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
    hsml = ngb_treefind(pos, hsml * 1.1);

    // save hsml for this search position
    hsml_out[i] = hsml;
  }

  tree_treefree();
  //free(P); // leave it alone

  IF_VERBOSE(printf("]\n"));
	
  return 1;
}

/* IDL glue routine
 
  example:

    nNGB = 32
    
    ; prepare inputs
    NumPart   = long(npts)
    Pos       = fltarr(3,npts)
    NumSrc    = long(nsearchpts)
    SearchPos = fltarr(3,nsearchpts)
    DesNumNgb    = long(nNGB)
    DesNumNgbDev = long(0)
    BoxSize      = 0.0
    
    hsml_out = fltarr(NumPart)
    
    ; call CalcHSMLds
    ret = Call_External('/n/home07/dnelson/idl/CalcHSMLds/CalcHSMLds_3D.so', 'CalcHSMLds', $
                        NumPart,Pos,NumSrc,SearchPos,DesNumNgb,DesNumNgbDev,BoxSize,hsml_out,/CDECL)
*/

int CalcHSMLds(int argc, void* argv[])
{
  int NumSearch;
  float *Pos,*SrcPos,*hsml_out;
	
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
  Pos       = (float *)argv[1];
  NumSearch = *(int *)argv[2];
  SrcPos    = (float *)argv[3];
	
  DesNumNgb    = *(int *)argv[4];
  DesNumNgbDev = *(int *)argv[5];
  BoxSize      = *(float *)argv[6];
  Softening    = 1.0; // code units
	
  hsml_out = (float *)argv[7];

#ifdef VERBOSE
  //printf("Input data:\n");
  //printf("NumPart      = %d\n", NumPart);
  //printf("NumSearch    = %d\n", NumSearch);
  //printf("DesNumNgb    = %d\n", DesNumNgb);
  //printf("DesNumNgbDev = %d\n", DesNumNgbDev);
  //printf("BoxSize      = %g\n", BoxSize);
#endif

  // calculate HSMLs
  return CalcHSMLds_natural(NumPart,Pos,NumSearch,SrcPos,DesNumNgb,DesNumNgbDev,BoxSize,hsml_out);
}
