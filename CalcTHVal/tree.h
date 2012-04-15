#ifndef TREE_H
#define TREE_H

#include "main.h"

#ifdef   ONEDIMS
#define  NUMDIMS 1
#define  KERNEL_COEFF_1  (4.0/3)
#define  KERNEL_COEFF_2  (8.0)
#define  KERNEL_COEFF_3  (24.0)
#define  KERNEL_COEFF_4  (16.0)
#define  KERNEL_COEFF_5  (8.0/3)
#define  KERNEL_COEFF_6  (-8.0)
#define  NORM_COEFF      2.0
#else
#ifndef  TWODIMS
#define  NUMDIMS 3    /**< For 3D-normalized kernel */
#define  KERNEL_COEFF_1  2.546479089470 /**< Coefficients for SPH spline kernel and its derivative */
#define  KERNEL_COEFF_2  15.278874536822
#define  KERNEL_COEFF_3  45.836623610466
#define  KERNEL_COEFF_4  30.557749073644
#define  KERNEL_COEFF_5  5.092958178941
#define  KERNEL_COEFF_6  (-15.278874536822)
#define  NORM_COEFF      4.188790204786 /**< Coefficient for kernel normalization. Note:  4.0/3 * PI = 4.188790204786 */
#else
#define  NUMDIMS 2    /**< For 2D-normalized kernel */
#define  KERNEL_COEFF_1  (5.0/7*2.546479089470) /**< Coefficients for SPH spline kernel and its derivative */
#define  KERNEL_COEFF_2  (5.0/7*15.278874536822)
#define  KERNEL_COEFF_3  (5.0/7*45.836623610466)
#define  KERNEL_COEFF_4  (5.0/7*30.557749073644)
#define  KERNEL_COEFF_5  (5.0/7*5.092958178941)
#define  KERNEL_COEFF_6  (5.0/7*(-15.278874536822))
#define  NORM_COEFF      M_PI /**< Coefficient for kernel normalization. */
#endif
#endif /* ONEDIMS */

/*! Macro that maps a distance to the nearest periodic neighbour */
#define NEAREST(x,boxsize) (((x)>boxhalf)?((x)-boxsize):(((x)<-boxhalf)?((x)+boxsize):(x)))
#define DMAX(a,b) (dmax1=(a),dmax2=(b),(dmax1>dmax2)?dmax1:dmax2)

struct NODE
{
  MyFloat len;
  MyFloat center[3];
  MyFloat maxsoft;
  union
  {
    int suns[8];
    struct
    {
      MyFloat s[3];
      MyFloat mass;
      unsigned int bitflags;
      int sibling;
      int nextnode;
      int father;
    }
    d;
  }
  u;
}
 *Nodes, *Nodes_base;

int MaxNodes;
int Numnodestree;
int TotNumNonInternalTopLevelTreeNodes;

int *Nextnode;
int *Prevnode;
int *Father;


#endif

