#ifndef __CSGRAPH_H__
#define __CSGRAPH_H__

#include <vector>

/*
 * Determine connected compoments of a compressed sparse graph.
 * Note:
 *   Output array flag must be preallocated
 */
template <class I>
I cs_graph_components(const I n_nod,
		      const I Ap[],
		      const I Aj[],
		            I flag[])
{
  // pos is a work array: list of nodes (rows) to process.
  std::vector<I> pos(n_nod,01);
  I n_comp = 0;
  I n_tot, n_pos, n_pos_new, n_pos0, n_new, n_stop;
  I icomp, ii, ir, ic;

  n_stop = n_nod;
  for (ir = 0; ir < n_nod; ir++) {
    flag[ir] = -1;
    if (Ap[ir+1] == Ap[ir]) {
      n_stop--;
      flag[ir] = -2;
    }
  }

  n_tot = 0;
  for (icomp = 0; icomp < n_nod; icomp++) {
    // Find seed.
    ii = 0;
    while ((flag[ii] >= 0) || (flag[ii] == -2)) {
      ii++;
      if (ii >= n_nod) {
	/* Sanity check, if this happens, the graph is corrupted. */
	return -1;
      }
    }

    flag[ii] = icomp;
    pos[0] = ii;
    n_pos0 = 0;
    n_pos_new = n_pos = 1;

    for (ii = 0; ii < n_nod; ii++) {
      n_new = 0;
      for (ir = n_pos0; ir < n_pos; ir++) {
	for (ic = Ap[pos[ir]]; ic < Ap[pos[ir]+1]; ic++) {
	  if (flag[Aj[ic]] == -1) {
	    flag[Aj[ic]] = icomp;
	    pos[n_pos_new] = Aj[ic];
	    n_pos_new++;
	    n_new++;
	  }
	}
      }
      n_pos0 = n_pos;
      n_pos = n_pos_new;
      if (n_new == 0) break;
    }
    n_tot += n_pos;

    if (n_tot == n_stop) {
      n_comp = icomp + 1;
      break;
    }
  }

  return n_comp;
}

#endif
#ifndef __CSGRAPH_H__
#define __CSGRAPH_H__

#include <vector>

/*
 * Determine connected compoments of a compressed sparse graph.
 * Note:
 *   Output array flag must be preallocated
 */
template <class I>
I cs_graph_components(const I n_nod,
		      const I Ap[],
		      const I Aj[],
		            I flag[])
{
  // pos is a work array: list of nodes (rows) to process.
  std::vector<I> pos(n_nod,01);
  I n_comp;
  I n_tot, n_pos, n_pos_new, n_pos0, n_new, n_stop;
  I icomp, ii, ir, ic;

  n_stop = n_nod;
  for (ir = 0; ir < n_nod; ir++) {
    flag[ir] = -1;
    if ((Ap[ir+1] - Ap[ir]) == 0) n_stop--;
  }

  n_tot = 0;
  for (icomp = 0; icomp < n_nod; icomp++) {
    // Find seed.
    ii = 0;
    while (flag[ii] >= 0) {
      ii++;
      if (ii >= n_nod) {
	/* Sanity check, if this happens, the graph is corrupted. */
	return -1;
      }
    }
    flag[ii] = icomp;
    pos[0] = ii;
    n_pos0 = 0;
    n_pos_new = n_pos = 1;

    for (ii = 0; ii < n_nod; ii++) {
      n_new = 0;
      for (ir = n_pos0; ir < n_pos; ir++) {
	for (ic = Ap[pos[ir]]; ic < Ap[pos[ir]+1]; ic++) {
	  if (flag[Aj[ic]] == -1) {
	    flag[Aj[ic]] = icomp;
	    pos[n_pos_new] = Aj[ic];
	    n_pos_new++;
	    n_new++;
	  }
	}
      }
      n_pos0 = n_pos;
      n_pos = n_pos_new;
      if (n_new == 0) break;
    }
    n_tot += n_pos;

    if (n_tot == n_stop) {
      n_comp = icomp + 1;
      break;
    }
  }

  return n_comp;
}

#endif
