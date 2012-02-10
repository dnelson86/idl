//
// CalcSphMap (sph kernel density estimation for projected maps) Routine for IDL
// dnelson jan.2012
// based on SphMap (for Python) by Mark Vogelsberger
//

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "idl_export.h"
#include "main.h"

/* project spline kernel */
inline double _getkernel(float h, float r2)
{
  double hinv, u;

  hinv = 1.0 / h;
  u = sqrt(r2) * hinv;

  if(u < 0.5)
    return (KERNEL_COEFF_1 + KERNEL_COEFF_2 * (u - 1.0) * u * u);
  else
    return KERNEL_COEFF_3 * (1.0 - u) * (1.0 - u) * (1.0 - u);

}

// main routine in C
int CalcSphMap_natural(float* data_pos, float* data_hsml, float* data_mass, float* grid_out)
{
  // working
  double pos0, pos1, h, v, /*w,*/ r2, h2;
  double p[3];
  double x, y;
  double pixelsizeX, pixelsizeY;
  double sum, pixelarea;
  double hmin, hmax;
  int nx, ny, dx, dy, i, j, k;
  double xx, yy, xxx, yyy;
  unsigned int binned_particles=0;
	
	// setup pixel sizes
  pixelsizeX = BoxSize_Lx / Xpixels;
  pixelsizeY = BoxSize_Ly / Ypixels;
  pixelarea = pixelsizeX * pixelsizeY;

  if(pixelsizeX < pixelsizeY)
    hmin = 1.001 * pixelsizeX / 2;
  else
    hmin = 1.001 * pixelsizeY / 2;

  hmax = pixelsizeX * 50;

#ifdef VERBOSE
  printf("calculating particle contributions...\n [");
  int signal = 0;
#endif
	
	// loop over all particles
  for (k=0; k < NumPart; k++)
  {
#ifdef VERBOSE
	  // output progress marker
		if (k > (signal / 100.0) * NumPart) {
			printf("x");
			fflush(stdout);
			signal++;
		}
#endif

    // get particle data
		p[0] = data_pos[3*k+0];
		p[1] = data_pos[3*k+1];
		p[2] = data_pos[3*k+2];
		h    = data_hsml[k];
		v    = data_mass[k];
		//w    = data_quantity[k];
		
		// check position, set (axis1,axis2) position, check hsml
		if(fabs(p[3 - axis0 - axis1] - BoxCen_Z) > BoxSize_Lz / 2.0)
			continue;

		pos0 = p[axis0] - (BoxCen_X - BoxSize_Lx / 2.0);
		pos1 = p[axis1] - (BoxCen_Y - BoxSize_Ly / 2.0);

		if(h < hmin)
			h = hmin;

		if(h > hmax)
			h = hmax;

		if(pos0 - 0 < -h || pos1 - 0 < -h || pos0 - BoxSize_Lx > h || pos1 - BoxSize_Ly > h)
			continue;

		binned_particles++;

		h2 = h * h;

		/* number of pixels covered by particle */
		nx = h / pixelsizeX + 1;
		ny = h / pixelsizeY + 1;

		/* coordinates of pixel center of particle */
		x = (floor(pos0 / pixelsizeX) + 0.5) * pixelsizeX;
		y = (floor(pos1 / pixelsizeY) + 0.5) * pixelsizeY;

		// calculate sum (normalization)
		sum = 0.0;

		for(dx = -nx; dx <= nx; dx++)
		{
			for(dy = -ny; dy <= ny; dy++)
			{
				/* distance of covered pixel from actual position */
				xx = x + dx * pixelsizeX - pos0;
				yy = y + dy * pixelsizeY - pos1;
				r2 = xx * xx + yy * yy;

				if(r2 < h2)
					sum += _getkernel(h, r2);
			}
		}

		if (sum < 1.0e-10)
			continue;

		// calculate contribution
		for(dx = -nx; dx <= nx; dx++)
		{
			for(dy = -ny; dy <= ny; dy++)
			{
				/* coordinates of pixel center of covering pixels */
				if(Periodic)
				{
					xxx = NEAREST(x + dx * pixelsizeX, BoxSize_Lx);
					yyy = NEAREST(y + dy * pixelsizeY, BoxSize_Ly);
				}
				else
				{
					xxx = x + dx * pixelsizeX;
					yyy = y + dy * pixelsizeY;
				}

				/* pixel array indices */
				i = xxx / pixelsizeX;
				j = yyy / pixelsizeY;

				if(i >= 0 && i < Xpixels)
				{
					if(j >= 0 && j < Ypixels)
					{
						xx = x + dx * pixelsizeX - pos0;
						yy = y + dy * pixelsizeY - pos1;
						r2 = xx * xx + yy * yy;

						if(r2 < h2)
						{
							// divide by sum for normalization
							grid_out[i * Ypixels + j] += _getkernel(h, r2) * v / sum;
							//gridquantity[i * Ypixels + j] += _getkernel(h, r2) * v * w / sum;
						}
					} //j
				} //i
			} //dy
		} //dx
		
  } //k

  IF_VERBOSE(printf("]\ndone.\n"));	
	
	// if requested (mode=3), divide by pixelarea to get column density map
  if(mode == 1) {
      return 1;
  }
  if(mode == 2) {
		printf("ERROR: QUANTITY MAP NOT IMPLEMENTED!\n");
		return 0;
  }

  if(mode == 3) {
		for(i = 0; i < Xpixels; i++)
			for(j = 0; j < Ypixels; j++)
				grid_out[i * Ypixels + j] /= pixelarea;

      return 1;
    }

  return 0;
}

// IDL glue routine
int CalcSphMap(int argc, void* argv[])
{	
  float *data_pos, *data_hsml, *data_mass, *grid_out;
	char buf[128];

  // validate input
  if (argc != 17)
  {
	  sprintf(buf,"Wrong number of arguments (%d)!\n",argc);
    IDL_Message(IDL_M_GENERIC,IDL_MSG_RET,buf);
    return 0;
  } else {
#ifdef TWODIMS
	  IDL_Message(IDL_M_GENERIC,IDL_MSG_RET,"CalcSphMap Loaded (TWODIMS!).");
#else
		IDL_Message(IDL_M_GENERIC,IDL_MSG_RET,"CalcSphMap Loaded (NDIMS=3!).");
#endif
  }
	
  // inputs and return by reference
  NumPart    = *(int *)argv[0];
  data_pos   = (float *)argv[1];
  data_hsml  = (float *)argv[2];
  data_mass  = (float *)argv[3];
	
	grid_out   = (float *)argv[4];
	
	// set globals
	BoxSize_Lx = *(float *)argv[5];
	BoxSize_Ly = *(float *)argv[6];
	BoxSize_Lz = *(float *)argv[7];
	
	BoxCen_X = *(float *)argv[8];
	BoxCen_Y = *(float *)argv[9];
	BoxCen_Z = *(float *)argv[10];
	
	axis0      = *(int *)argv[11];
	axis1      = *(int *)argv[12];
	
	Xpixels    = *(int *)argv[13];
	Ypixels    = *(int *)argv[14];
	mode       = *(int *)argv[15];
	Periodic   = *(int *)argv[16];

#ifdef VERBOSE
  printf("Input data:\n");
	printf("NumPart     = %d\n", NumPart);
  printf("BoxSize     = %g %g %g\n", BoxSize_Lx, BoxSize_Ly, BoxSize_Lz);
  printf("BoxCen      = %g %g %g\n", BoxCen_X, BoxCen_Y, BoxCen_Z);
  printf("axes        = %d %d\n", axis0, axis1);
  printf("nPixels     = %d %d\n", Xpixels, Ypixels);
  printf("mode per    = %d %d\n", mode, Periodic);
#endif

	// calculate density projection map
	return CalcSphMap_natural(data_pos,data_hsml,data_mass,grid_out);
}