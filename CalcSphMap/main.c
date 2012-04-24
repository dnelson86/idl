//
// CalcSphMap (sph kernel density estimation for projected maps) Routine for IDL
// dnelson apr.2012
// based on SphMap (for Python) by Mark Vogelsberger (based on HsmlAndProject for IDL)
// note: changed to column major output to match IDL convention
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
int CalcSphMap_natural(float* data_pos, float* data_hsml, float* data_mass, float* data_quant,
                       float* dens_out, float* quant_out)
{
  // working
  double pos0, pos1, h, v, w, r2, h2;
  double p[3];
  double x, y;
  double pixelsizeX, pixelsizeY;
  double sum, pixelarea;
  double hmin, hmax;
  int nx, ny, dx, dy, i, j, k;
  double xx, yy, xxx, yyy;
  unsigned int binned_particles=0;

  // setup pixel sizes
  pixelsizeX = BoxSize_ImageX / Xpixels;
  pixelsizeY = BoxSize_ImageY / Ypixels;
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
    w    = data_quant[k];
		
    // clip points outside box (z) dimension
    if(BoxSize)
    {
      if(fabs(NEAREST(p[3-axis0-axis1]-BoxCen_Z)) > 0.5*BoxSize_ImageZ)
        continue;
    } else {
      if(fabs(p[3 - axis0 - axis1] - BoxCen_Z) > 0.5*BoxSize_ImageZ)
        continue;
    }

    // position relative to box (x,y) minimum
    pos0 = p[axis0] - (BoxCen_X - 0.5*BoxSize_ImageX);
    pos1 = p[axis1] - (BoxCen_Y - 0.5*BoxSize_ImageY);

    // clamp smoothing length
    if(h < hmin)
      h = hmin;

    if(h > hmax)
      h = hmax;

    // clip points outside box (x,y) dimensions
    if(BoxSize) // periodic
    {
      if(NEAREST(BoxCen_X-p[axis0]) + 0.5*BoxSize_ImageX < -h || 
         NEAREST(BoxCen_Y-p[axis1]) + 0.5*BoxSize_ImageY < -h || 
         NEAREST(BoxCen_X-p[axis0]) - 0.5*BoxSize_ImageX > +h || 
         NEAREST(BoxCen_Y-p[axis1]) - 0.5*BoxSize_ImageY > +h)
        continue;
    } else {
      if(pos0 - 0 < -h || pos1 - 0 < -h || pos0 - BoxSize_ImageX > h || pos1 - BoxSize_ImageY > h)
        continue;
    }

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
        if(BoxSize)
	{
	  xxx = NEAREST(x + dx * pixelsizeX);
	  yyy = NEAREST(y + dy * pixelsizeY);
          //xxx = NEAREST(xxx - 0.5*BoxSize_ImageX,BoxSize) + 0.5*BoxSize_ImageX;
          //yyy = NEAREST(yyy - 0.5*BoxSize_ImageY,BoxSize) + 0.5*BoxSize_ImageY;
	} else {
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
		dens_out[j * Ypixels + i]  += _getkernel(h, r2) * v / sum;
		quant_out[j * Ypixels + i] += _getkernel(h, r2) * v * w / sum;
            }
          } //j
        } //i
      } //dy
    } //dx
		
  } //k

  // normalize mass weighted quantity
  for(j = 0; j < Ypixels; j++)
    for(i = 0; i < Xpixels; i++)
      if(dens_out[j * Ypixels + i] > 0)
        quant_out[j * Ypixels + i] /= dens_out[j * Ypixels + i];

  IF_VERBOSE(printf("]\n"));	
	
  return 1;
}

// IDL glue routine
int CalcSphMap(int argc, void* argv[])
{	
  float *data_pos, *data_hsml, *data_mass, *data_quant;
  float *dens_out, *quant_out;
  char buf[128];

  // validate input
  if (argc != 18)
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
  data_quant = (float *)argv[4];
	
  dens_out   = (float *)argv[5];
  quant_out  = (float *)argv[6];
  
  // set globals
  BoxSize_ImageX = *(float *)argv[7];
  BoxSize_ImageY = *(float *)argv[8];
  BoxSize_ImageZ = *(float *)argv[9];
  
  BoxSize    = *(float *)argv[10];
  BoxHalf    = BoxSize / 2.0;
  BoxCen_X   = *(float *)argv[11];
  BoxCen_Y   = *(float *)argv[12];
  BoxCen_Z   = *(float *)argv[13];
	
  axis0      = *(int *)argv[14];
  axis1      = *(int *)argv[15];
	
  Xpixels    = *(int *)argv[16];
  Ypixels    = *(int *)argv[17];

#ifdef VERBOSE
  //printf("Input data:\n");
  //printf("NumPart     = %d\n", NumPart);
  //printf("BoxSize_Im  = %g %g %g\n", BoxSize_ImageX, BoxSize_ImageY, BoxSize_ImageZ);
  //printf("BoxSize     = %g\n", BoxSize);
  //printf("BoxCen      = %g %g %g\n", BoxCen_X, BoxCen_Y, BoxCen_Z);
  //printf("axes        = %d %d\n", axis0, axis1);
  //printf("nPixels     = %d %d\n", Xpixels, Ypixels);
#endif

  // calculate density projection map
  return CalcSphMap_natural(data_pos,data_hsml,data_mass,data_quant,dens_out,quant_out);
}
