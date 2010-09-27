
0.4
===

Major changes in this release include:

   - Coordinate Descent algorithm (Lasso, ElasticNet) refactoring & 
     speed improvements (roughly 100x times faster).

   - Coordinate Descent Refactoring (and bug fixing) for consistency
     with R's package GLMNET.

   - New metrics module.

   - New GMM module contributed by Ron Weiss.

   - Implementation of the LARS algorithm (without Lasso variant for now).

   - feature_selection module redesign.

   - Migration to GIT as content management system.

   - Removal of obsolete attrselect module.

   - Rename of private compiled extensions (aded underscore).

   - Removal of legacy unmaintained code.

   - Documentation improvements (both docstring and rst).

   - Improvement of the build system to (optionally) link with MKL. 
Also, provide a lite BLAS implementation in case no system-wide BLAS is 
found.

   - Lots of new examples.

   - Many, many bug fixes ...


Authors
-------

The committer list for this release is the following (preceded by number 
of commits):

    * 143  Fabian Pedregosa
    * 35  Alexandre Gramfort
    * 34  Olivier Grisel
    * 11  Gael Varoquaux
    *  5  Yaroslav Halchenko
    *  2  Vincent Michel
    *  1  Chris Filo Gorgolewski

