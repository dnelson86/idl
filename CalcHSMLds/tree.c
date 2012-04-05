//
// CalcHSMLds Routine for IDL
// dnelson apr.2012
// tree.c based on CalcHsml (for Python) by Mark Vogelsberger
//

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#define  MAX_REAL_NUMBER  1.0e35

//#include "main.h"
#include "tree.h"

static int last;		/* auxialiary variable used to set-up non-recursive walk */

int
tree_treebuild (void)
{
  int i, j, subnode = 0, parent = -1, numnodes;
  int nfree, no, nn;
  double lenhalf;
  struct NODE *nfreep;
  float xmin[3], xmax[3], len;


  Nodes = Nodes_base - NumPart;

  /* select first node */
  nfree = NumPart;
  nfreep = &Nodes[nfree];

  /* create an empty root node  */
  if (BoxSize>0.0) /* PERIODIC */
    {
      for (j = 0; j < 3; j++)
	nfreep->center[j] = 0.5 * BoxSize;
      nfreep->len = BoxSize;
    }
  else /* NON PERIODIC */
    {
      for (j = 0; j < 3; j++)
	{
	  xmin[j] = MAX_REAL_NUMBER;
	  xmax[j] = -MAX_REAL_NUMBER;
	}

      for (i = 0; i < NumPart; i++)
	{
	  for (j = 0; j < 3; j++)
	    {
	      if (P[i].Pos[j] > xmax[j])
		xmax[j] = P[i].Pos[j];
	      if (P[i].Pos[j] < xmin[j])
		xmin[j] = P[i].Pos[j];
	    }
	}

      /* determine maxmimum extension */
      len = xmax[0] - xmin[0];
      if ((xmax[1] - xmin[1]) > len)
	len = xmax[1] - xmin[1];
      if ((xmax[2] - xmin[2]) > len)
	len = xmax[2] - xmin[2];

      for (j = 0; j < 3; j++)
	nfreep->center[j] = 0.5 * (xmin[j] + xmax[j]);
      nfreep->len = len;
    }

  /* daugther slots are all empty */
  for (i = 0; i < 8; i++)
    nfreep->u.suns[i] = -1;

  numnodes = 1;
  nfreep++;
  nfree++;
 
  /* at this point the root node with empty daugthers is created */

  /* insert all particles / contruct tree */
  for (i = 0; i < NumPart; i++)
    {
      /* start at the root node */
      no = NumPart;

      /* insert particle i */
      while (1)
	{
	  if (no >= NumPart)	/* we are dealing with an internal node */
	    {
	      /* to which subnode will this particle belong */
	      subnode = 0;
	      if (P[i].Pos[0] > Nodes[no].center[0])
		subnode += 1;
	      if (P[i].Pos[1] > Nodes[no].center[1])
		subnode += 2;
	      if (P[i].Pos[2] > Nodes[no].center[2])
		subnode += 4;

              /* get the next node */ 
	      nn = Nodes[no].u.suns[subnode];

	      if (nn >= 0)	/* ok, something is in the daughter slot already, need to continue */
		{
		  parent = no;	/* note: subnode can still be used in the next step of the walk */
		  no = nn;
		}
	      else
		{
		  /* here we have found an empty slot where we can 
		   * attach the new particle as a leaf 
		   */
		  Nodes[no].u.suns[subnode] = i;
		  break;	/* done for this particle */
		}
	    }
	  else
	    {
	      /* we try to insert into a leaf with a single particle
	       * need to generate a new internal node at this point 
	       because every leaf is only allowed to contain one particle
	       */
	      Nodes[parent].u.suns[subnode] = nfree;

	      nfreep->len = 0.5 * Nodes[parent].len;
	      lenhalf = 0.25 * Nodes[parent].len;


	      if (subnode & 1)
		nfreep->center[0] = Nodes[parent].center[0] + lenhalf;
	      else
		nfreep->center[0] = Nodes[parent].center[0] - lenhalf;

	      if (subnode & 2)
		nfreep->center[1] = Nodes[parent].center[1] + lenhalf;
	      else
		nfreep->center[1] = Nodes[parent].center[1] - lenhalf;

	      if (subnode & 4)
		nfreep->center[2] = Nodes[parent].center[2] + lenhalf;
	      else
		nfreep->center[2] = Nodes[parent].center[2] - lenhalf;

	      nfreep->u.suns[0] = -1;
	      nfreep->u.suns[1] = -1;
	      nfreep->u.suns[2] = -1;
	      nfreep->u.suns[3] = -1;
	      nfreep->u.suns[4] = -1;
	      nfreep->u.suns[5] = -1;
	      nfreep->u.suns[6] = -1;
	      nfreep->u.suns[7] = -1;

	      subnode = 0;
	      if (P[no].Pos[0] > nfreep->center[0])
		subnode += 1;
	      if (P[no].Pos[1] > nfreep->center[1])
		subnode += 2;
	      if (P[no].Pos[2] > nfreep->center[2])
		subnode += 4;

	      if (nfreep->len < 1.0e-3 * Softening)
		{
		  /* seems like we're dealing with particles at identical locations.   
		   * Randomize subnode index to get particle in different leaf.
		   * (happens well below gravitational softening scale). 
		   */
		  subnode = (int) (8.0 * drand48 ());
		  if (subnode >= 8)
		    subnode = 7;
		}

	      nfreep->u.suns[subnode] = no;

	      no = nfree;	/* resume trying to insert the new particle at
				   the newly created internal node */

	      numnodes++;
	      nfree++;
	      nfreep++;

	      if ((numnodes) >= MaxNodes)
		{
		  printf ("maximum number %d of tree-nodes reached.\n",
			  MaxNodes);
		  printf ("for particle %d  %g %g %g\n", i, P[i].Pos[0],
			  P[i].Pos[1], P[i].Pos[2]);
		  exit (1);
		}
	    }
	}
    }

  /* now compute the multipole moments recursively */

  last = -1;
  tree_update_node_recursive (NumPart, -1, -1);

  if (last >= NumPart)
    Nodes[last].u.d.nextnode = -1;
  else
    Nextnode[last] = -1;


  return numnodes;
}

/* this routine computes the multipole moments for a given internal node and
 * all its subnodes using a recursive computation.  Note that the moments of
 * the daughter nodes are already stored in single precision. For very large
 * particle numbers, loss of precision may results for certain particle
 * distributions
 */
void
tree_update_node_recursive (int no, int sib, int father)
{
  int j, jj, p, pp = 0, nextsib, mass, suns[8];
  double s[3];

  if (no >= NumPart)
    {
      for (j = 0; j < 8; j++)
	suns[j] = Nodes[no].u.suns[j];	/* this "backup" is necessary because the nextnode entry will
					   overwrite one element (union!) */
      if (last >= 0)
	{
	  if (last >= NumPart)
	    Nodes[last].u.d.nextnode = no;
	  else
	    Nextnode[last] = no;
	}

      last = no;

      mass = 0;
      s[0] = 0;
      s[1] = 0;
      s[2] = 0;

      for (j = 0; j < 8; j++)
	{
	  if ((p = suns[j]) >= 0)
	    {
	      /* check if we have a sibling on the same level */
	      for (jj = j + 1; jj < 8; jj++)
		if ((pp = suns[jj]) >= 0)
		  break;

	      if (jj < 8)	/* yes, we do */
		nextsib = pp;
	      else
		nextsib = sib;

	      tree_update_node_recursive (p, nextsib, no);

	      if (p >= NumPart)	/* an internal node or pseudo particle */
		{
		  mass += Nodes[p].u.d.mass;	/* we assume a fixed particle mass */
		  s[0] += Nodes[p].u.d.mass * Nodes[p].u.d.s[0];
		  s[1] += Nodes[p].u.d.mass * Nodes[p].u.d.s[1];
		  s[2] += Nodes[p].u.d.mass * Nodes[p].u.d.s[2];
		}
	      else		/* a particle */
		{
		  mass += 1.;
		  s[0] += P[p].Pos[0];
		  s[1] += P[p].Pos[1];
		  s[2] += P[p].Pos[2];
		}
	    }
	}

      if (mass)
	{
	  s[0] /= mass;
	  s[1] /= mass;
	  s[2] /= mass;
	}
      else
	{
	  s[0] = Nodes[no].center[0];
	  s[1] = Nodes[no].center[1];
	  s[2] = Nodes[no].center[2];
	}

      Nodes[no].u.d.s[0] = s[0];
      Nodes[no].u.d.s[1] = s[1];
      Nodes[no].u.d.s[2] = s[2];
      Nodes[no].u.d.mass = mass;

      Nodes[no].u.d.sibling = sib;
      Nodes[no].u.d.father = father;
    }
  else				/* single particle or pseudo particle */
    {
      if (last >= 0)
	{
	  if (last >= NumPart)
	    Nodes[last].u.d.nextnode = no;
	  else
	    Nextnode[last] = no;
	}

      last = no;

      if (no < NumPart)		/* only set it for single particles */
	Father[no] = father;
    }
}

float ngb_treefind (float xyz[3], float hguess)
{
  int iter;
  float numngb;
  float left=0, right=0;
  float dmax1, dmax2;

  if (hguess == 0)
    {
      hguess = 1.;
    }

  iter = 0;

  do
    {
      iter++;
      if (iter>1000)
       {
        printf("too many iterations...\n");
        exit(1);
       }

      numngb = ngb_treefind_variable (xyz, hguess);
      if (numngb < (DesNumNgb-DesNumNgbDev) || numngb >= (DesNumNgb+DesNumNgbDev))
      {


      if (left > 0 && right > 0)
       if (right-left < 0.001 * left)
        break; // particle is OK

      if(numngb < (DesNumNgb - DesNumNgbDev))
       left = DMAX(hguess, left);
      else
       {
        if(right != 0)
         {
          if(hguess < right)
           right = hguess;
          }
        else
         right = hguess;
       }

                  if(right > 0 && left > 0)
                    hguess = pow(0.5 * (pow(left, 3) + pow(right, 3)), 1.0 / 3);
                  else
                    {
                      if(right == 0 && left == 0)
                        exit(1);   /* can't occur */

                      if(right == 0 && left > 0)
                        hguess *= 1.26;

                      if(right > 0 && left == 0)
                        hguess /= 1.26;
                    }

      }
     else
      break;

    }
  while (1);
  return hguess;
}

float ngb_treefind_variable (float searchcenter[3], float hguess)
{
  int numngb, no, p;
  double dx, dy, dz, r2, h2, r, weighted_numngb, wk, hinv, hinv3, u;
  struct NODE *this;

  double boxhalf = 0.5 * BoxSize;

  h2 = hguess * hguess;
  hinv = 1.0/hguess;
  hinv3 = hinv*hinv*hinv;
  
  numngb = 0;
  weighted_numngb = 0.0;
  no = NumPart;

  while (no >= 0)
    {
      if (no < NumPart)		/* single particle */
	{
	  p = no;
	  no = Nextnode[no];

	  dx = NEAREST (P[p].Pos[0] - searchcenter[0], BoxSize);
	  if (dx < -hguess)
	    continue;
	  if (dx > hguess)
	    continue;

	  dy = NEAREST (P[p].Pos[1] - searchcenter[1], BoxSize);
	  if (dy < -hguess)
	    continue;
	  if (dy > hguess)
	    continue;

	  dz = NEAREST (P[p].Pos[2] - searchcenter[2], BoxSize);
	  if (dz < -hguess)
	    continue;
	  if (dz > hguess)
	    continue;

	  r2 = dx * dx + dy * dy + dz * dz;

	  if (r2 < h2)
	    {
                  r = sqrt(r2);

                  u = r * hinv;

                  if(u < 0.5)
                    {
                      wk = hinv3 * (KERNEL_COEFF_1 + KERNEL_COEFF_2 * (u - 1) * u * u);
                    }
                  else
                    {
                      wk = hinv3 * KERNEL_COEFF_5 * (1.0 - u) * (1.0 - u) * (1.0 - u);
                    } 


                  weighted_numngb += NORM_COEFF * wk / hinv3;     /* 4.0/3 * PI = 4.188790204786 */
	      numngb++;
	    }
	}
      else
	{
	  this = &Nodes[no];

	  no = Nodes[no].u.d.sibling;	/* in case the node can be discarded */

	  if ((NEAREST (this->center[0] - searchcenter[0], BoxSize) +
	       0.5 * this->len) < -hguess)
	    continue;
	  if ((NEAREST (this->center[0] - searchcenter[0], BoxSize) -
	       0.5 * this->len) > hguess)
	    continue;
	  if ((NEAREST (this->center[1] - searchcenter[1], BoxSize) +
	       0.5 * this->len) < -hguess)
	    continue;
	  if ((NEAREST (this->center[1] - searchcenter[1], BoxSize) -
	       0.5 * this->len) > hguess)
	    continue;
	  if ((NEAREST (this->center[2] - searchcenter[2], BoxSize) +
	       0.5 * this->len) < -hguess)
	    continue;
	  if ((NEAREST (this->center[2] - searchcenter[2], BoxSize) -
	       0.5 * this->len) > hguess)
	    continue;

	  no = this->u.d.nextnode;	/* ok, we need to open the node */
	}
    }


  return weighted_numngb;
}


/* this function allocates memory used for storage of the tree
 * and auxiliary arrays for tree-walk and link-lists.
 */
size_t
tree_treeallocate (int maxnodes, int maxpart)	/* usually maxnodes=0.7*maxpart is sufficient */
{
  size_t bytes, allbytes = 0;

  MaxNodes = maxnodes;

  if (!(Nodes_base = malloc (bytes = (MaxNodes + 1) * sizeof (struct NODE))))
    {
      printf ("failed to allocate memory for %d tree-nodes (%g MB).\n",
	      MaxNodes, bytes / (1024.0 * 1024.0));
      exit (3);
    }
  allbytes += bytes;

  if (!(Nextnode = malloc (bytes = maxpart * sizeof (int))))
    {
      printf ("Failed to allocate %d spaces for 'Nextnode' array (%g MB)\n",
	      maxpart, bytes / (1024.0 * 1024.0));
      exit (4);
    }
  allbytes += bytes;

  if (!(Father = malloc (bytes = maxpart * sizeof (int))))
    {
      printf ("Failed to allocate %d spaces for 'Father' array (%g MB)\n",
	      maxpart, bytes / (1024.0 * 1024.0));
      exit (5);
    }
  allbytes += bytes;

  return allbytes;
}


/* free the allocated memory
 */
void
tree_treefree (void)
{
  free (Father);
  free (Nextnode);
  free (Nodes_base);
}
