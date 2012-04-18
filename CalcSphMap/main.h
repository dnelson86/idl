#ifndef MAIN_H
#define MAIN_H

#define VERBOSE

//#define NEAREST(x) (((x)>BoxSize)?((x)-BoxSize):(((x)<0)?((x)+BoxSize):(x))) //no
#define NEAREST(x) (((x)>BoxHalf)?((x)-BoxSize):(((x)<-BoxHalf)?((x)+BoxSize):(x))) //yes

#ifdef VERBOSE
#define IF_VERBOSE(expr) (expr)
#else
#define IF_VERBOSE(expr) ((void)0)
#endif

int NumPart;
double BoxSize;
double BoxHalf;
float BoxSize_ImageX, BoxSize_ImageY, BoxSize_ImageZ;
float BoxCen_X, BoxCen_Y, BoxCen_Z;
int axis0, axis1;
int Xpixels, Ypixels;

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

