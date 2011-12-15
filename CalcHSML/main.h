#ifndef ALLVARS_H
#define ALLVARS_H

//#define VERBOSE

//#define TWODIMS // affects kernel coefficients, set in makefile
//#define STRIDE_BYTES_PER_ELEM 4 // 4 byte floats, unused

#ifdef VERBOSE
#define IF_VERBOSE(expr) (expr)
#else
#define IF_VERBOSE(expr) ((void)0)
#endif

// tree.c definitions
int tree_treebuild(void);
float ngb_treefind(float[],float);
float ngb_treefind_variable(float[],float);
size_t tree_treeallocate(int,int);
void tree_treefree(void);
void tree_update_node_recursive(int,int,int);

// globals
int NumPart;
int DesNumNgb;
int DesNumNgbDev;
float BoxSize;
float Softening;

// single precision
typedef unsigned int MyIDType;
typedef float MyFloat;
typedef float MyDouble;

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


#endif
