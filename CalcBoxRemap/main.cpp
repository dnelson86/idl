//
// CalcBoxRemap
// dnelson oct.2013
// "A volume and local structure preserving mapping of periodic boxes"
// based on: http://mwhite.berkeley.edu/BoxRemap/
//

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "idl_export.h"
#include "cuboid.h"

/* IDL glue routine
 
  example:

  ; prepare inputs
	pos         = fltarr(3,npos)
  NumPart     = long64(npos)
  BoxSize     = float(boxsize)
  remapMatrix = long(remapMatrix) ; 3x3 integer matrix, row-major

  ret = Call_External('/n/home07/dnelson/idl/CalcBoxRemap/CalcBoxRemap.so', 'CalcBoxRemap', $
                      NumPart,Pos,BoxSize,remapMatrix,/CDECL)
*/

extern "C" int CalcBoxRemap(int argc, void* argv[]) /* note extern "C" to disable c++ name mangling */
{
  long long NumPart;
	double BoxSize;
	float *Pos;
	int remapMatrix[9];
	char buf[128];
	int signal = 0;

  // validate input
  if (argc != 4)
  {
	  sprintf(buf,"Wrong number of arguments (%d)!\n",argc);
    IDL_Message(IDL_M_GENERIC,IDL_MSG_RET,buf);
    return 0;
  }

  // inputs and return by reference
  NumPart = *(long long *)argv[0];
  Pos     = (float *)argv[1];
	BoxSize = *(float *)argv[2];
	
	int *matElement = (int *)argv[3];
	
	for(int i = 0; i < 9; i++)
	  remapMatrix[i] = *(matElement + i);

#ifdef VERBOSE
  printf("Input data:\n");
	printf("NumPart      = %lld\n", NumPart);
  printf("BoxSize      = %g\n", BoxSize);
	for(int i=0; i < 9; i++)
    printf("Matrix[%d]     = %d\n", i, remapMatrix[i]);
#endif

  // initialize mapping
	Cuboid R(remapMatrix);
	
	// intermediate variables: do computations in double
	double x1, x2, x3;
  double r1, r2, r3;
	double invBoxSize = 1.0 / BoxSize;
	
	printf("boxRemap: [");
	
	for(long long i=0; i < NumPart; i++)
	{
#ifdef VERBOSE
		if(i < 5)
			printf("[%lld] In:  %g %g %g\n", i, Pos[ 3*i + 0 ], Pos[ 3*i + 1 ], Pos[ 3*i + 2 ]);
#endif
		if (i > (signal / 100.0) * NumPart)
		{
			printf("x");
			fflush(stdout);
			signal++;
		}
			
		// convert input points to [0,1]^3 domain (in double)
		x1 = Pos[ 3*i + 0 ] * invBoxSize;
		x2 = Pos[ 3*i + 1 ] * invBoxSize;
		x3 = Pos[ 3*i + 2 ] * invBoxSize;

		// transform points into new cuboid
		R.Transform(x1, x2, x3, r1, r2, r3);
		
		// convert back to [0,L1],[0,L2],[0,L3] domain (in float)
#ifdef SKIPZ_STRIDE2
		Pos[ 2*i + 0 ] = r1 * BoxSize;
		Pos[ 2*i + 1 ] = r2 * BoxSize;
#else
		Pos[ 3*i + 0 ] = r1 * BoxSize;
		Pos[ 3*i + 1 ] = r2 * BoxSize;
		Pos[ 3*i + 2 ] = r3 * BoxSize;
#endif

#ifdef VERBOSE
		if(i < 5)
			printf("[%lld] Out: %g %g %g\n", i, Pos[ 3*i + 0 ], Pos[ 3*i + 1 ], Pos[ 3*i + 2 ]);
#endif
	}
	
	printf("]\n");
	
	return 1;
}

/* Inverse (single point only, input/output normalized in [0,L1/2/3] or [0,1]
 
  example:

  ; prepare inputs
	px = double(px)
	py = double(py)
	pz = double(pz)
  remapMatrix = long(remapMatrix) ; 3x3 integer matrix, row-major

  ret = Call_External('/n/home07/dnelson/idl/CalcBoxRemap/CalcBoxRemap.so', 'CalcBoxRemapInv', $
                      px,py,pz,remapMatrix,/CDECL)
*/

extern "C" int CalcBoxRemapInv(int argc, void* argv[]) /* note extern "C" to disable c++ name mangling */
{
	double *r1,*r2,*r3;
	int remapMatrix[9];
	char buf[128];

  // validate input
  if (argc != 4)
  {
	  sprintf(buf,"Wrong number of arguments (%d)!\n",argc);
    IDL_Message(IDL_M_GENERIC,IDL_MSG_RET,buf);
    return 0;
  }

  // inputs and return by reference
  r1 = (double *)argv[0];
  r2 = (double *)argv[1];
	r3 = (double *)argv[2];
	
	int *matElement = (int *)argv[3];
	
	for(int i = 0; i < 9; i++)
	  remapMatrix[i] = *(matElement + i);

#ifdef VERBOSE
  printf("Input data:\n");
  printf("r1      = %g\n", *r1);
  printf("r2      = %g\n", *r2);
  printf("r3      = %g\n", *r3);
	for(int i=0; i < 9; i++)
    printf("Matrix[%d]     = %d\n", i, remapMatrix[i]);
#endif

  // initialize mapping
	Cuboid R(remapMatrix);
	
  double x1, x2, x3;
	
	// transform points into new cuboid
	R.InverseTransform(*r1, *r2, *r3, x1, x2, x3);
	
	// overwrite inputs
	*r1 = x1;
	*r2 = x2;
	*r3 = x3;
	
	return 1;
}
