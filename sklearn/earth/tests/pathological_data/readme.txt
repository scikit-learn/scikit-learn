Pathological Data Sets

The data sets contained in this folder have revealed bugs in the past.  They now serve as 
regression tests to prevent the same bugs from arising again in the future.

issue_44:
This data set caused a segfault during fitting due to the use of a negative index in Cython
with wraparound = False.