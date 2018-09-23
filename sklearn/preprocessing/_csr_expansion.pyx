# cython: cdivision=True
# cython: boundscheck=False
# cython: wraparound=False

# Author: Andrew nystrom <awnystrom@gmail.com>

# When testing, test the following cases:
# A single dense row
# An empty matrix
# Multiple rows, but one is all zeros.
# The first row is all zeros. Tests we're making indptr correctly.

cimport numpy as np
cimport cython
from scipy.sparse import csr_matrix
from numpy cimport ndarray

ctypedef np.float64_t DATA_T
ctypedef np.int32_t INDEX_T
ctypedef np.int64_t INDEX_T64

def csr_expansion_deg2(ndarray[DATA_T, ndim=1] data not None,
         ndarray[INDEX_T, ndim=1] indices not None,
         ndarray[INDEX_T, ndim=1] indptr not None,
         INDEX_T D, INDEX_T interaction_only):
    """
    Perform a second-degree polynomial or interaction expansion on a compressed scipy
    sparse row (CSR) matrix.
    
    Parameters
    ----------
    data : the "data" attribute of the input CSR matrix.
    indices : the "indices" attribute of the input CSR matrix.
    indptr : the "indptr" attribute of the input CSR matrix.
    D : The dimensionality of the input CSR matrix.
    interaction_only : 0 for a polynomial expansion, 1 for an interaction expansion.
    """
    
    cdef INDEX_T total_nnz = 0, row_i, nnz_in_row
  
    cdef INDEX_T expanded_this_row
    # Count how many nonzero elements the expanded matrix will contain.
    for row_i in range(indptr.shape[0]-1):
      nnz_in_row = indptr[row_i + 1] - indptr[row_i]
      if interaction_only:
        expanded_this_row = <INDEX_T>((nnz_in_row**2 - nnz_in_row) / 2)
      else:
        expanded_this_row = <INDEX_T>((nnz_in_row**2 + nnz_in_row) / 2)
      total_nnz += expanded_this_row
  
    # Make the arrays that will form the CSR matrix of the expansion.
    cdef ndarray[DATA_T, ndim=1] expanded_data = ndarray(shape=total_nnz, dtype='float64',
                                                         order='C')
    cdef ndarray[INDEX_T, ndim=1] expanded_indices = ndarray(shape=total_nnz,
                                                               dtype='int32', order='C')
    cdef ndarray[INDEX_T, ndim=1] expanded_indptr = ndarray(shape=indptr.shape[0],
                                                              dtype='int32', order='C')
  
    expanded_indptr[0] = indptr[0]
    cdef INDEX_T expanded_index = 0, row_starts, row_ends, i, j, k1, k2, \
      num_cols_in_row, expanded_column_part
    
    #Calculate the augmented matrix and polynomial expansion features.
    for row_i in range(indptr.shape[0]-1):
      row_starts = indptr[row_i]
      row_ends = indptr[row_i + 1]
      num_cols_in_row = 0
      for k1 in range(row_starts, row_ends):
        i = indices[k1]
        # Cache this part as it doesn't involve j.
        if interaction_only:
          expanded_column_part = D*i - (i**2 + 3*i) / 2 - 1
        else:
          expanded_column_part = D*i - (i**2 + i) / 2
        #Add the expansion features.
        for k2 in range(k1 + interaction_only, row_ends):
          j = indices[k2]
          #Now put everything in its right place.
          expanded_data[expanded_index] = data[k1] * data[k2]
          expanded_indices[expanded_index] = expanded_column_part + j
          expanded_index += 1
          num_cols_in_row += 1
    
      expanded_indptr[row_i+1] = expanded_indptr[row_i] + num_cols_in_row
      
    X_poly = csr_matrix([])
    X_poly.data = expanded_data
    X_poly.indices = expanded_indices
    X_poly.indptr = expanded_indptr
    if interaction_only:
      X_poly._shape = (indptr.shape[0]-1, <INDEX_T>((D**2 - D) / 2))
    else:
      X_poly._shape = (indptr.shape[0]-1, <INDEX_T>((D**2 + D) / 2))
    return X_poly