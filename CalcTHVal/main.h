#ifndef ALLVARS_H
#define ALLVARS_H

#define VERBOSE

#ifdef VERBOSE
#define IF_VERBOSE(expr) (expr)
#else
#define IF_VERBOSE(expr) ((void)0)
#endif

#define PI 3.1415926535897932384626434
#define PI_INV 0.31830988618379067153

// tree.c definitions
int tree_treebuild(void);
float ngb_treefind(float[],float,float*);
float ngb_treefind_variable(float[],float,int*);
size_t tree_treeallocate(int,int);
void tree_treefree(void);
void tree_update_node_recursive(int,int,int);

// globals
int NumPart;
int DesNumNgb;
int DesNumNgbDev;
int TophatMode;
float BoxSize;
float Softening;

int *Ngblist;

// single precision
//typedef unsigned int MyIDType;
typedef float MyFloat;
typedef float MyDouble;

// particle_data subset from allvars.h (minimal)
struct particle_data
{
  MyDouble Pos[3];		/**< particle position at its current time */
  MyFloat Value;                /* value of some function specified for each particle */
  //MyDouble Vel[3];		/**< particle velocity at its current time */
  //MyDouble Mass;		/**< particle mass */
	
  //MyIDType ID;
  //short int Type;		/**< flags particle type.  0=gas, 1=halo, 2=disk, 3=bulge, 4=stars, 5=bndry */
}
 *P;


#endif
