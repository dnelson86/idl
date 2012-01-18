#ifndef MAIN_H
#define MAIN_H

#define VERBOSE

//#define TWODIMS // affects kernel coefficients, set in makefile

#define NEAREST(x, BoxSize) (((x)>BoxSize)?((x)-BoxSize):(((x)<0)?((x)+BoxSize):(x)))

#ifdef VERBOSE
#define IF_VERBOSE(expr) (expr)
#else
#define IF_VERBOSE(expr) ((void)0)
#endif

int NumPart;
float BoxSize_Lx, BoxSize_Ly, BoxSize_Lz;
float BoxCen_X, BoxCen_Y, BoxCen_Z;
int axis0, axis1;
int Xpixels, Ypixels, mode, Periodic;

#ifdef   ONEDIMS
#define  NUMDIMS 1
#define  KERNEL_COEFF_1  (4.0/3)
#define  KERNEL_COEFF_2  (8.0)
#define  KERNEL_COEFF_3  (24.0)
#else
#ifndef  TWODIMS
#define  NUMDIMS 3    /**< For 3D-normalized kernel */
#define  KERNEL_COEFF_1  2.546479089470 /**< Coefficients for SPH spline kernel and its derivative */
#define  KERNEL_COEFF_2  15.278874536822
#define  KERNEL_COEFF_3  45.836623610466
#else
#define  NUMDIMS 2    /**< For 2D-normalized kernel */
#define  KERNEL_COEFF_1  (5.0/7*2.546479089470) /**< Coefficients for SPH spline kernel and its derivative */
#define  KERNEL_COEFF_2  (5.0/7*15.278874536822)
#define  KERNEL_COEFF_3  (5.0/7*45.836623610466)
#endif
#endif /* ONEDIMS */

#endif

