Support for Python 3
====================

There is experimental support for python3 in some modules, the status
of these is:

   svm : OK
   linear_models : ?
   cluster : ?
   covariance : ?
   grid_search : ?
   mixture : ?
   externals.joblib : ?
   utils.sparsetools : FAILS to compile

To generate python3 compatible sources for selected modules, run the
2to3 tool on the module::

    2to3 -wn --no-diffs scikits/learn/$module

If you would like to help with porting to python3, please propose
yourself in the scikit-learn mailing list:

    https://lists.sourceforge.net/lists/listinfo/scikit-learn-general
