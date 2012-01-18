//
// CalcHSML Routine for IDL
// dnelson dec.2011
// tree.c taken from CalcHsml (for Python) by Mark Vogelsberger
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

/** This function returns the gas particle nearest to the specified
    position by searching the tree. It only looks among the local
    particles, and only looks within the specified mindist. If no gas
    particle is within the set mindist, -1 is returned. Upon return,
    mindist is set to the actual distance to the nearest particle. The
    tree traversal starts at startnode. To start at the root node, set
    it to -1.

    The initial guess for mindist has a severe impact on the search
    efficiency, so a good guess is desired. If no guess is known, set
    mindist and guess to -1 and a random gas particle will be picked as the
    initial guess. The initial guess can be specified by either a gas
    particle in guess or a distance in mindist.

    If the search is to use periodic boundary conditions, set boxSize > 0.
 */
int ngb_treefind_nearest_local(MyDouble searchcenter[3], int startnode, MyDouble * mindist, int guess)
{
  int node, nearest, p, niter = 0;
  struct NODE *current;
  MyDouble dx, dy, dz;
  MyFloat cur_mindist = *mindist;

#ifdef PERIODIC
  MyDouble xtmp;
#endif
  if(startnode < 0)
    node = NumPart;
  else
    node = startnode;

  /** it's better to have a guess of the min distance, because this
   minimizes the number of node openings we have to do. If it was not
   provided, let's just pick a random gas particle. */
  if(guess >= 0)
    {
      // override mindist with guess
      *mindist = -1;
      nearest = guess;
    }
  else
    //nearest = floor(get_random_number(SelRnd++) * N_gas);
    nearest = floor(NumPart/2);

  if(*mindist < 0)
    {
      if(boxSize>0.0)
	{
	  dx = NGB_PERIODIC_LONG_X(P[nearest].Pos[0] - searchcenter[0]);
	  dy = NGB_PERIODIC_LONG_Y(P[nearest].Pos[1] - searchcenter[1]);
	  dz = NGB_PERIODIC_LONG_Z(P[nearest].Pos[2] - searchcenter[2]);
	}
      else
	{
	  dx = fabs(P[nearest].Pos[0] - searchcenter[0]);
	  dy = fabs(P[nearest].Pos[1] - searchcenter[1]);
	  dz = fabs(P[nearest].Pos[2] - searchcenter[2]);
	}
      cur_mindist = sqrt(dx * dx + dy * dy + dz * dz);
    }
  else
    {
      nearest = -1;
      cur_mindist = *mindist;
    }

  int endnode = Nodes[node].u.d.sibling;
  while(node != endnode)
    {
      ++niter;

      if(node < NumPart)	/* single particle */
	{
	  p = node;
	  node = Nextnode[node];

	  if(P[p].Type > 0)
	    {
	      /* Not gas particle.  This can happen because there are
	         several particles hanging off the smallest nodes, and we
	         only know that one of them must be a gas particle from
	         our flag test. */
	      continue;
	    }

	  if(boxSize>0.0)
	    dx = NGB_PERIODIC_LONG_X(P[p].Pos[0] - searchcenter[0]);
	  else
	    dx = fabs(P[p].Pos[0] - searchcenter[0]);
	  if(dx > cur_mindist)
	    continue;
	  if(boxSize>0.0)
	    dy = NGB_PERIODIC_LONG_Y(P[p].Pos[1] - searchcenter[1]);
	  else
	    dy = fabs(P[p].Pos[1] - searchcenter[1]);
	  if(dy > cur_mindist)
	    continue;
	  if(boxSize>0.0)
	    dz = NGB_PERIODIC_LONG_Z(P[p].Pos[2] - searchcenter[2]);
	  else
	    dz = fabs(P[p].Pos[2] - searchcenter[2]);
	  if(dz > cur_mindist)
	    continue;
	  MyDouble curdist2 = dx * dx + dy * dy + dz * dz;
	  if(curdist2 > cur_mindist * cur_mindist)
	    continue;

	  cur_mindist = sqrt(curdist2);
	  nearest = p;
	}
      else
	{
	  if(node >= NumPart + MaxNodes)
	    {
	      // hitting a pseudo particle means the search box
	      // extended into another task domain and we can no
	      // longer guarantee that the nearest particle is on this
	      // processor.
	      // SHOULD NOT HAPPEN HERE
        printf("ERROR1\n");
	      node = Nextnode[node - MaxNodes];
	      continue;
	    }

	  current = &Nodes[node];

	  node = current->u.d.sibling;	/* in case the node can be discarded */

	  // first quick tests along the axes
	  MyDouble test_dist = cur_mindist + 0.5 * current->len;

	  if(boxSize>0.0)
	    dx = NGB_PERIODIC_LONG_X(current->center[0] - searchcenter[0]);
	  else
	    dx = fabs(current->center[0] - searchcenter[0]);
	  if(dx > test_dist)
	    continue;

	  if(boxSize>0.0)
	    dy = NGB_PERIODIC_LONG_Y(current->center[1] - searchcenter[1]);
	  else
	    dy = fabs(current->center[1] - searchcenter[1]);
	  if(dy > test_dist)
	    continue;

	  if(boxSize>0.0)
	    dz = NGB_PERIODIC_LONG_Z(current->center[2] - searchcenter[2]);
	  else
	    dz = fabs(current->center[2] - searchcenter[2]);
	  if(dz > test_dist)
	    continue;

	  /* now test against the minimal sphere enclosing everything */
	  test_dist += FACT1 * current->len;
	  if(dx * dx + dy * dy + dz * dz > test_dist * test_dist)
	    continue;

	  node = current->u.d.nextnode;	/* ok, we need to open the node */
	}
    }
#ifdef VERBOSE
  //printf("nearest particle %d found (dist=%f) in %d iterations\n", nearest, *mindist, niter);
#endif
  *mindist = cur_mindist;
  return nearest;
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
