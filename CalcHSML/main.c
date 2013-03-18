//
// CalcHSML Routine for IDL
// dnelson feb.2012
// based on CalcHsml (for Python) by Mark Vogelsberger
// 1D,2D,3D versions, for periodic set BoxSize>0 (LONG_XYZ not supported though)
//

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "idl_export.h"
#include "main.h"

// main routine in C
int CalcHSML_natural(int NumPart, float* data_pos, float* data_mass, int DesNumNgb, int DesNumngbDev, 
                     float BoxSize, float HsmlGuess, float Softening, float* hsml_out)
{
  int i;
  int nNodes;
  float hsml = HsmlGuess;

  P = (struct particle_data *) malloc (NumPart * sizeof (struct particle_data));
	
  // copy data into P
  for (i=0; i < NumPart; i++)
  {
	P[i].Pos[0] = data_pos[3*i+0];
	P[i].Pos[1] = data_pos[3*i+1];
	P[i].Pos[2] = data_pos[3*i+2];
	P[i].Mass   = data_mass[i];
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
	
  for (i=0; i < NumPart; i++)
  {
#ifdef VERBOSE
	  // output progress marker
		if (i > (signal / 100.0) * NumPart)
		{
			printf("x");
			fflush(stdout);
			signal++;
		}
#endif
		// calculate smoothing length and save
		hsml = ngb_treefind(P[i].Pos, hsml * 1.1);
		
		hsml_out[i] = hsml;
  }

  tree_treefree();
  free(P);

  IF_VERBOSE(printf("]\ndone.\n"));
	
  return 1;
}

/* IDL glue routine
 
  example:

    nNGB = 32
    
    ; prepare inputs
    NumPart = long(npts)
    Pos     = fltarr(3,npts)
    Mass    = fltarr(npts)
    
    DesNumNgb    = long(nNGB)
    DesNumNgbDev = long(0)
    BoxSize      = 0.0
    HsmlGuess    = float(1.0)
    Softening    = float(1.0)
    
    hsml_out = fltarr(NumPart)
    
    ; call CalcHSML
    ret = Call_External('/n/home07/dnelson/idl/CalcHSML/CalcHSML.so', 'CalcHSML', $
                        NumPart,Pos,Mass,DesNumNgb,DesNumNgbDev,BoxSize,HsmlGuess,Softening,hsml_out, $
                        /CDECL)
*/

int CalcHSML(int argc, void* argv[])
{
	float HsmlGuess;
  float *Pos,*Mass,*hsml_out;
	
	char buf[128];

  // validate input
  if (argc != 9)
  {
	  sprintf(buf,"Wrong number of arguments (%d)!\n",argc);
    IDL_Message(IDL_M_GENERIC,IDL_MSG_RET,buf);
    return 0;
  } else {
#ifdef VERBOSE
#ifdef   ONEDIMS
IDL_Message(IDL_M_GENERIC,IDL_MSG_RET,"CalcHSML Loaded (ONEDIMS!).");
#else
#ifndef  TWODIMS
IDL_Message(IDL_M_GENERIC,IDL_MSG_RET,"CalcHSML Loaded (NDIMS=3!).");
#else
IDL_Message(IDL_M_GENERIC,IDL_MSG_RET,"CalcHSML Loaded (TWODIMS!).");
#endif
#endif /* ONEDIMS */
#endif
  }

  // should check: pos has 2 dims, first dim is 3
  // should check: mass has 1 dim
	
  // inputs and return by reference
  NumPart = *(int *)argv[0];
  Pos     = (float *)argv[1];
  Mass    = (float *)argv[2];
	
  DesNumNgb    = *(int *)argv[3];
	DesNumNgbDev = *(int *)argv[4];
	BoxSize      = *(float *)argv[5];
	HsmlGuess    = *(float *)argv[6];
	Softening    = *(float *)argv[7];
	
	hsml_out = (float *)argv[8];

#ifdef VERBOSE
  printf("Input data:\n");
	printf("NumPart      = %d\n", NumPart);
  printf("DesNumNgb    = %d\n", DesNumNgb);
  printf("DesNumNgbDev = %d\n", DesNumNgbDev);
  printf("BoxSize      = %g\n", BoxSize);
  printf("HsmlGuess    = %g\n", HsmlGuess);
  printf("Softening    = %g\n", Softening);
#endif

	// calculate HSMLs
	return CalcHSML_natural(NumPart,Pos,Mass,DesNumNgb,DesNumNgbDev,BoxSize,HsmlGuess,Softening,hsml_out);
}
