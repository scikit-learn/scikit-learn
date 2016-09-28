# Ryan Peach 3/1/2016
# References Used for this Implementation
# https://en.wikipedia.org/wiki/Hungarian_algorithm
# https://github.com/bmc/munkres/blob/master/munkres.py
# http://csclab.murraystate.edu/bob.pilgrim/445/munkres.html
# ---
# No copying and pasting of online code was performed, though some code may turn out to be similar due to the standardized nature of this algorithm.
# This is a very different implementation due to the fact that it heavily uses numpy, vastly simplifies many poorly pythonized coding elements
# removes the "class" approach for a function based one. Etc.
# Improvements that need to be made is to require the matrix class, and to vectorize iterations through the matrix.

import numpy as np
from queue import PriorityQueue
from time import time
import functools

# Module Globals
NONE, ZERO, STAR, PRIME = -1, 0, 1, 2
hungmem = {}

def score(C,answ):
    return sum([C[n,answ[n]] for n in range(len(answ))])

def hungarian(C):
    """ Calculates the hungarian of C, remembering all values to avoid recalculation. """
    global hungmem            # Use module memory
    #print(C)
    C = np.matrix(C)          # Use a matrix
    H = C.tostring()
    if H not in hungmem:      # Check if it's in memory
        S = _hungarian(C)     # If it's not, calculate the solution
        hungmem[H] = S        # Then add the solution to memory
    else:                     # If it is
        S = hungmem[H]        # simply retrieve the solution from memory
    return S                  # Return the answer

def murty(P0):
    """ Non-optimized Murty's. Generator.
        Ref: Optimizing Murty's Ranked Assignment Method, Fig. 4
        by Matt L. Miller, Harold S. Stone, & Ingemar J. Cox """
    try:
        INF = np.iinfo(P0.dtype).max                   # A value used to remove (y,z,l) from a problem
    except:
        INF = np.finfo(P0.dtype).max

    def valid(X):
        """ Checks if any rows or columns are completely ignored. """
        for row in np.array(X):    # Iterate over all rows
            found = False
            for y in row:          # Iterate over all values in row
                if y != INF:       # If you find a cell that is not ignored
                    found = True
                    break          # Go on to the next row
            if not found:
                return False       # If you never find one, this is an invalid matrix

        for col in np.array(X.T):  # Iterate over all columns
            found = False
            for x in col:          # Iterate over all values in column
                if x != INF:       # If you find a cell that is not ignored
                    found = True
                    break          # Go on to the next column
            if not found:
                return False       # If you never find one, this is an invalid matrix

        return True

    S0 = hungarian(P0)                              # Find the best solution S0 to problem P0 (1)
    C0 = score(P0, S0)                              # Find the cost of S0 given problem P0
    Q  = PriorityQueue()                            # Initialize a priority queue (C, P, S) sorted top as lowest cost (2)
    Q.put((C0,time(),P0.copy(),S0))                 # Add (C0, P0, S0) to the queue. Use time to sort equal costs
                                                    # P0 will be modified, but the original is needed for score references.
    solutions, found = set(), set()                 # Create a searchable set of solutions. Found used to avoid duplicates (3)
    while not Q.empty():                            # Iterate until all solutions found (4)
        #print("Murtys")
        C, _, P, S = Q.get()                        # Get the top set off the Queue (4.1)
        solutions.add(tuple(S))                     # Add solution to output set for searching (4.2)
        yield (S,score(P0,S))                       # The generator for this function returns each solution in order.
        for y, x in enumerate(S):                   # Iterating over the solution, ignoring the last (4.3)
            l = P0[y,x]                             # l is the cost of the assignment y to z
            newP = P.copy()                         # copy P (4.3.1)
            newP[y,x] = INF                         # Remove y, z from newP (4.3.2)
            if valid(newP):                         # if new solution exists (4.3.4)
                newS = hungarian(newP)              # find new solution to new problem (4.3.3)
                if tuple(newS) not in found:        # avoid duplicate answers
                    newC = score(P0,newS)           # get the cost of newS in P0
                    Q.put((newC,time(),newP,newS))  # Add (newC, newP, newS) to the queue, use time to sort equal costs
                    found.add(tuple(newS))          # Add newP, newS to found
            for i, j in np.ndindex(*P.shape):       # From P, remove y row and z column, except for index y, z
                if (i == x and j != y) or (i != x and j == y):
                    P[j,i] = INF

    raise StopIteration()

def _hungarian(C):
    """ Calculates the hungarian of the cost matrix C.
                This function performs steps 0, 1, and 2 of the munkres method.
        Complexity: O(N**3) where N is the dimensionality of C
        Time: 0.52N**3 - 2.57N**2 + 104.63N - 111.74 (ms) (See Hungarian.ods)
        Ref: http://csclab.murraystate.edu/bob.pilgrim/445/munkres.html """
    return _Hstep0(np.matrix(C).copy().T)  # calls step 0 on a copy of matrix C, transversed to make output correlate with assigned columns


def _Hstep0(C):
    """ Create matrix such that #columns >= #rows """
    if not C.shape[1] >= C.shape[0]:  # If there are not more columns than rows, traspose the matrix
        C = C.T
    return _Hstep1(C)

def _Hstep1(C):
    """ For each row, find the minimum and subtract it from it's row. """
    C -= C.min(axis=1)
    return _Hstep2(C)

def _Hstep2(C):
    """ For all zeroes, star zero if no other stared zero's exist in it's row or col. """
    global NONE, ZERO, STAR, PRIME
    marked = np.matrix(np.zeros(C.shape))
    ycovered, xcovered = [False for y in range(C.shape[0])], [False for x in range(C.shape[1])]
    for y in range(C.shape[0]):
        for x in range(C.shape[1]):
            if (C[y,x] == ZERO) and not xcovered[x] and not ycovered[y]:
                marked[y, x] = STAR
                xcovered[x] = True
                ycovered[y] = True
    ycovered, xcovered = [False for y in range(C.shape[0])], [False for x in range(C.shape[1])]

    return _Hstep3(C,marked,xcovered,ycovered)

def _Hstep3(C,marked,xcovered,ycovered):
    """ Step 3:  Creates covers. Returns final answer. """
    global NONE, ZERO, STAR, PRIME
    for y in range(C.shape[0]):
        for x in range(C.shape[1]):
            if marked[y,x] == STAR and not xcovered[x] and not ycovered[y]:  # Cover columns which contain a star
                xcovered[x] = True
    if all(xcovered):                                                        # If all columns are covered, we have found our answer
        return [findM(col,STAR)[1] for col in marked.T]                      # Returun the row indexes of all STARS for every column
    else:
        #print("Step 3: ",C,ycovered,xcovered,marked)
        return _Hstep4(C,marked,xcovered,ycovered)                             # Otherwise, proceed to step 4

def _Hstep4(C,marked,xcovered,ycovered):
    """ Step 4:  Changes the covers to cover Primes and uncover STARS.
                Primes uncovered zeros. """
    global NONE, ZERO, STAR, PRIME
    while True:
        (y,x) = findM(C,ZERO,xcovered,ycovered)                 # Find first ZERO in covered Cost Matrix
        if y != NONE:
            marked[y,x] = PRIME                                 # Prime the uncovered zero
            if STAR not in marked[y,:]:                         # If there is no STAR in the row
                #print("Step 4: ",C,ycovered,xcovered,marked)
                return _Hstep5(C,y,x,marked,xcovered,ycovered)  # Go to step 5
            else:                                               # Otherwise
                ycovered[y] = True                              # Cover the row
                xcovered[findM(marked[y,:],STAR)[1]] = False    # Uncover the column containing the found STAR
        else:
            #print("Step 4: ",C,ycovered,xcovered,marked)
            return _Hstep6(C,ycovered,xcovered,marked)            # Go to step 6 if there is no ZERO in Cost Matrix

def _Hstep5(C,y,x,marked,xcovered,ycovered):
    """ Step 5:  Finds a path stairstepping from first star to first prime in star's row to first prime in star's column, etc.
                Similar to Stable Marriage Algorithm. """
    global NONE, ZERO, STAR, PRIME
    path = []
    path.append((y,x))
    while True:
        r = findM(marked[:,path[-1][1]],STAR)[0]        # Find row with a star
        if r != NONE:
            path.append((r,path[-1][1]))                # Append row found alongside last column found
        else:
            break                                       # Done when none found
        c = findM(marked[path[-1][0],:],PRIME)[1]       # Find column with a prime
        path.append((path[-1][0],c))                    # Append column found alongside last row found

    # Unstar each starred zero of the series, star each primed zero of the series.
    for y, x in path:
        if marked[y,x] == STAR:
            marked[y,x] = ZERO
        else:
            marked[y,x] = STAR

    # Erase all primes
    erasure = np.vectorize(lambda x: ZERO if x == PRIME else x)
    marked = erasure(marked)

    # Uncover every line in the matrix
    ycovered, xcovered = [False for y in range(C.shape[0])], [False for x in range(C.shape[1])]

    # Return to step 3
    #print("Step 5: ",C,ycovered,xcovered,marked,path)
    return _Hstep3(C,marked,xcovered,ycovered)

def _Hstep6(C,ycovered,xcovered,marked):
    """ Step 6:  Add the min value to each covered rows, and subtract it from all uncovered columns. """
    # Calculate the minimum uncovered value
    minv = uncovered(C,ycovered,xcovered)[0].min()

    # Add the min value to each covered rows, and subtract it from all uncovered columns.
    add6 = np.vectorize(lambda v,y: v+minv if ycovered[y] else v)
    sub6 = np.vectorize(lambda v,x: v-minv if not xcovered[x] else v)
    X, Y = np.meshgrid(np.arange(C.shape[1]),np.arange(C.shape[0]))
    C = add6(C,Y)
    C = sub6(C,X)

    #print("Step 6: ",C,ycovered,xcovered,marked)
    return _Hstep4(C,marked,xcovered,ycovered)

def uncovered(M, ycovered, xcovered):
    """ Returns a matrix identical to M but with the covered rows and columns deleted """
    M = M.copy()                                                      # Copy M so that it is not directly modified
    X, Y = np.meshgrid(np.arange(M.shape[1]),np.arange(M.shape[0]))
    xcovered = [x for x in range(len(xcovered)) if xcovered[x]]       # Enumerate indexes where x is covered
    ycovered = [y for y in range(len(ycovered)) if ycovered[y]]       # Enumerate indexes where y is covered
    M = np.delete(M, ycovered, axis=0)                                # Delete indexes where y is covered
    Y = np.delete(Y, ycovered, axis=0)
    M = np.delete(M, xcovered, axis=1)                                # Delete indexes where x is covered
    X = np.delete(X, xcovered, axis=1)
    return M, Y, X

def findM(M, val, xcovered = None, ycovered = None):
    """ Finds the first instance of val in M in rows and columns which are uncovered. """
    # Create default x and y covers
    if xcovered == None:
        xcovered = [False for x in range(M.shape[1])]
    if ycovered == None:
        ycovered = [False for y in range(M.shape[0])]

    # Find y and x where M[y,x] == val
    for y in range(M.shape[0]):
        for x in range(M.shape[1]):
            if not xcovered[x] and not ycovered[y]:
                if M[y,x] == val:
                    return y, x

    # Return global None if none found
    global NONE
    return NONE, NONE

# ------------ Test Methods ---------------
def timeHungarian(n2, n1 = 3):
    times = []
    for n in range(n1, n2):
        t = timeit('hungarian(np.random.rand({0},{0}))'.format(n), setup='from __main__ import hungarian; import numpy as np;')
        times.append(t)
    return times,range(n1,n2)

if __name__ == "__main__":
    # Test 1: hungarian([[1,2,3],[2,4,6],[3,6,9]]) == [2,1,0]
    A = np.matrix([[1,2,3],[2,4,6],[3,6,9]])
    S = hungarian(A)
    success = (S == [2,1,0])
    print("Test1 (Hungarian): " + str(success))
    print(str(A) + "\nHungarian -> " + str(S))

    # Test 2: Murty's Algorithm
    for S in murty(A):
        print(S)
