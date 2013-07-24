"""
Solve the unique lowest-cost assignment problem using the
Hungarian algorithm (also known as Munkres algorithm).

"""
# Based on original code by Brain Clapper, adapted to numpy by Gael Varoquaux

# Copyright (c) 2008 Brian M. Clapper <bmc@clapper.org>, Gael Varoquaux
# Author: Brian M. Clapper, Gael Varoquaux
# LICENSE: BSD

import numpy as np


###############################################################################
# Object-oriented form of the algorithm
class _Hungarian(object):
    """Hungarian algorithm

    Calculate the Munkres solution to the classical assignment problem.

    Warning: this code is not following scikit-learn standards and will be
    refactored.
    """

    def compute(self, cost_matrix):
        """
        Compute the indices for the lowest-cost pairings.

        Parameters
        ----------
        cost_matrix : 2D matrix
            The cost matrix. Does not have to be square.

        Returns
        -------
        indices : 2D array of indices
            The pairs of (row, col) indices in the original array giving
            the original ordering.
        """
        cost_matrix = np.atleast_2d(cost_matrix)

        # If there are more rows (n) than columns (m), then the algorithm
        # will not be able to work correctly. Therefore, we
        # transpose the cost function when needed. Just have to
        # remember to swap the result columns later in this function.
        doTranspose = (cost_matrix.shape[1] < cost_matrix.shape[0])
        if doTranspose:
            self.C = (cost_matrix.T).copy()
        else:
            self.C = cost_matrix.copy()

        # At this point, m >= n.
        self.n = n = self.C.shape[0]
        self.m = m = self.C.shape[1]
        self.row_uncovered = np.ones(n, dtype=np.bool)
        self.col_uncovered = np.ones(m, dtype=np.bool)
        self.Z0_r = 0
        self.Z0_c = 0
        self.path = np.zeros((n+m, 2), dtype=int)
        self.marked = np.zeros((n, m), dtype=int)

        done = False
        step = 1

        steps = {1: self._step1,
                 3: self._step3,
                 4: self._step4,
                 5: self._step5,
                 6: self._step6}

        if m == 0 or n == 0:
            # No need to bother with assignments if one of the dimensions
            # of the cost matrix is zero-length.
            done = True

        while not done:
            try:
                func = steps[step]
                step = func()
            except KeyError:
                done = True

        # Look for the starred columns
        results = np.array(np.where(self.marked == 1)).T

        # We need to swap the columns because we originally
        # did a transpose on the input cost matrix.
        if doTranspose:
            results = results[:, ::-1]

        return results.tolist()

    def _step1(self):
        """ Steps 1 and 2 in the wikipedia page.
        """
        # Step1: For each row of the matrix, find the smallest element and
        # subtract it from every element in its row.
        self.C -= self.C.min(axis=1)[:, np.newaxis]
        # Step2: Find a zero (Z) in the resulting matrix. If there is no
        # starred zero in its row or column, star Z. Repeat for each element
        # in the matrix.
        for i, j in zip(*np.where(self.C == 0)):
            if self.col_uncovered[j] and self.row_uncovered[i]:
                self.marked[i, j] = 1
                self.col_uncovered[j] = False
                self.row_uncovered[i] = False

        self._clear_covers()
        return 3

    def _step3(self):
        """
        Cover each column containing a starred zero. If n columns are
        covered, the starred zeros describe a complete set of unique
        assignments. In this case, Go to DONE, otherwise, Go to Step 4.
        """
        marked = (self.marked == 1)
        self.col_uncovered[np.any(marked, axis=0)] = False

        if marked.sum() >= self.n:
            return 7  # done
        else:
            return 4

    def _step4(self):
        """
        Find a noncovered zero and prime it. If there is no starred zero
        in the row containing this primed zero, Go to Step 5. Otherwise,
        cover this row and uncover the column containing the starred
        zero. Continue in this manner until there are no uncovered zeros
        left. Save the smallest uncovered value and Go to Step 6.
        """
        # We convert to int as numpy operations are faster on int
        C = (self.C == 0).astype(np.int)
        covered_C = C*self.row_uncovered[:, np.newaxis]
        covered_C *= self.col_uncovered.astype(np.int)
        n = self.n
        m = self.m
        while True:
            # Find an uncovered zero
            row, col = np.unravel_index(np.argmax(covered_C), (n, m))
            if covered_C[row, col] == 0:
                return 6
            else:
                self.marked[row, col] = 2
                # Find the first starred element in the row
                star_col = np.argmax(self.marked[row] == 1)
                if not self.marked[row, star_col] == 1:
                    # Could not find one
                    self.Z0_r = row
                    self.Z0_c = col
                    return 5
                else:
                    col = star_col
                    self.row_uncovered[row] = False
                    self.col_uncovered[col] = True
                    covered_C[:, col] = C[:, col] * (
                        self.row_uncovered.astype(np.int))
                    covered_C[row] = 0

    def _step5(self):
        """
        Construct a series of alternating primed and starred zeros as
        follows. Let Z0 represent the uncovered primed zero found in Step 4.
        Let Z1 denote the starred zero in the column of Z0 (if any).
        Let Z2 denote the primed zero in the row of Z1 (there will always
        be one). Continue until the series terminates at a primed zero
        that has no starred zero in its column. Unstar each starred zero
        of the series, star each primed zero of the series, erase all
        primes and uncover every line in the matrix. Return to Step 3
        """
        count = 0
        path = self.path
        path[count, 0] = self.Z0_r
        path[count, 1] = self.Z0_c
        done = False
        while not done:
            # Find the first starred element in the col defined by
            # the path.
            row = np.argmax(self.marked[:, path[count, 1]] == 1)
            if not self.marked[row, path[count, 1]] == 1:
                # Could not find one
                done = True
            else:
                count += 1
                path[count, 0] = row
                path[count, 1] = path[count-1, 1]

            if not done:
                # Find the first prime element in the row defined by the
                # first path step
                col = np.argmax(self.marked[path[count, 0]] == 2)
                if self.marked[row, col] != 2:
                    col = -1
                count += 1
                path[count, 0] = path[count-1, 0]
                path[count, 1] = col

        # Convert paths
        for i in range(count+1):
            if self.marked[path[i, 0], path[i, 1]] == 1:
                self.marked[path[i, 0], path[i, 1]] = 0
            else:
                self.marked[path[i, 0], path[i, 1]] = 1

        self._clear_covers()
        # Erase all prime markings
        self.marked[self.marked == 2] = 0
        return 3

    def _step6(self):
        """
        Add the value found in Step 4 to every element of each covered
        row, and subtract it from every element of each uncovered column.
        Return to Step 4 without altering any stars, primes, or covered
        lines.
        """
        # the smallest uncovered value in the matrix
        if np.any(self.row_uncovered) and np.any(self.col_uncovered):
            minval = np.min(self.C[self.row_uncovered], axis=0)
            minval = np.min(minval[self.col_uncovered])
            self.C[np.logical_not(self.row_uncovered)] += minval
            self.C[:, self.col_uncovered] -= minval
        return 4

    def _find_prime_in_row(self, row):
        """
        Find the first prime element in the specified row. Returns
        the column index, or -1 if no starred element was found.
        """
        col = np.argmax(self.marked[row] == 2)
        if self.marked[row, col] != 2:
            col = -1
        return col

    def _clear_covers(self):
        """Clear all covered matrix cells"""
        self.row_uncovered[:] = True
        self.col_uncovered[:] = True


###############################################################################
# Functional form for easier use
def linear_assignment(X):
    """Solve the linear assignment problem using the Hungarian algorithm

    The problem is also known as maximum weight matching in bipartite graphs.
    The method is also known as the Munkres or Kuhn-Munkres algorithm.

    Parameters
    ----------
    X : array
        The cost matrix of the bipartite graph

    Returns
    -------
    indices : array,
        The pairs of (row, col) indices in the original array giving
        the original ordering.

    References
    ----------

    1. http://www.public.iastate.edu/~ddoty/HungarianAlgorithm.html

    2. Harold W. Kuhn. The Hungarian Method for the assignment problem.
       *Naval Research Logistics Quarterly*, 2:83-97, 1955.

    3. Harold W. Kuhn. Variants of the Hungarian method for assignment
       problems. *Naval Research Logistics Quarterly*, 3: 253-258, 1956.

    4. Munkres, J. Algorithms for the Assignment and Transportation Problems.
       *Journal of the Society of Industrial and Applied Mathematics*,
       5(1):32-38, March, 1957.

    5. http://en.wikipedia.org/wiki/Hungarian_algorithm
    """
    H = _Hungarian()
    indices = H.compute(X)
    indices.sort()
    # Re-force dtype to ints in case of empty list
    indices = np.array(indices, dtype=int)
    # Make sure the array is 2D with 2 columns.
    # This is needed when dealing with an empty list
    indices.shape = (-1, 2)
    return indices
