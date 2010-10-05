
.. currentmodule:: scikits.learn

0.5
===

.. _changes_0_5:

Changelog
---------

New modules
~~~~~~~~~~~

    - New Hidden Markov Models module (:ref:`scikits.learn.hmm`)

    - Support for sparse matrices in modules ``svm`` and ``glm``.

    - Add a Pipeline object to compose different estimators.

    - Recursive Feature Elimination routines in module
      :ref:`modules/feature_selection`.

    - Addition of various classes capable of cross validation in the
      glm module (:class:`glm.LassoCV`, :class:`glm.ElasticNetCV`,
      etc.).

    - New pure-python LARS algorithm implementation
      (scikits.learn.glm.lars).


Examples
~~~~~~~~

    - new examples using some of the mlcomp datasets:
      :ref:`example_mlcomp_sparse_document_classification.py`,
      :ref:`example_mlcomp_document_classification.py`

    - Many more examaples.

Fixes
~~~~~

    - API changes: adhere variable names to PEP-8, give more
      meaningful names.



External dependencies
~~~~~~~~~~~~~~~~~~~~~

    - Joblib is now a dependencie of this package, although it is
      shipped with (scikits.learn.externals.joblib).

Removed modules
~~~~~~~~~~~~~~~

    - Module ann (Artificial Neural Networks) has been removed from
      the distribution. Users wanting this sort of algorithms should
      take a look into pybrain.

0.4
===

Changelog
---------

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

