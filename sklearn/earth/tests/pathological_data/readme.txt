Pathological Data Sets

The data sets contained in this folder have revealed bugs in the past.  They now serve as 
regression tests to prevent the same bugs from arising again in the future.

issue_44:
This data set caused a segfault during fitting due to the use of a negative index in Cython
with wraparound = False.

issue_50:
This data set exposed a bug that occurred when using the sample_weight parameter.  The problem 
was that the apply methods of the BasisFunctions were not applying the weights, and the 
next_pair method of the ForwardPasser class assumed they were.  Now next_pair applies the 
weights after calling apply.  The same data set exposed issue 51, in which user-specified 
endspans were not used.  This test case covers both issues.