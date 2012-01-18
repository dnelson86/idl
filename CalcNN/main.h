#ifndef ALLVARS_H
#define ALLVARS_H

#define VERBOSE

//#define TWODIMS // affects kernel coefficients, set in makefile

#define PERIODIC

#ifdef VERBOSE
#define IF_VERBOSE(expr) (expr)
#else
#define IF_VERBOSE(expr) ((void)0)
#endif

// single precision
typedef unsigned int MyIDType;
typedef float MyFloat;
typedef float MyDouble;

// tree.c definitions
int tree_treebuild(void);
size_t tree_treeallocate(int,int);
void tree_treefree(void);
void tree_update_node_recursive(int,int,int);

// based on newer ngb.c functions
int ngb_treefind_nearest_local(MyDouble[],int,MyDouble*,int);

// globals
int NumPart;
int DesNumNgb;
int DesNumNgbDev;
float BoxSize;
float Softening;

// particle_data subset from allvars.h
struct particle_data
{
  MyDouble Pos[3];		/**< particle position at its current time */
  MyDouble Vel[3];		/**< particle velocity at its current time */
  MyDouble Mass;		/**< particle mass */
	
  MyIDType ID;
  short int Type;		/**< flags particle type.  0=gas, 1=halo, 2=disk, 3=bulge, 4=stars, 5=bndry */
}
 *P;

#ifdef PERIODIC
//extern MyDouble boxSize, boxHalf;

#define boxSize BoxSize
#define boxHalf BoxSize*0.5

#define boxSize_X boxSize
#define boxHalf_X boxHalf
#define boxSize_Y boxSize
#define boxHalf_Y boxHalf
#define boxSize_Z boxSize
#define boxHalf_Z boxHalf
#endif

#endif
