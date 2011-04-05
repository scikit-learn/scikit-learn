
.. currentmodule:: scikits.learn


.. _changes_0_7:

0.7
===

scikits.learn 0.7 was released in March 2011, roughly three months
after the 0.6 release. This release is marked by the speed
improvements in existing algorithms like k-Nearest Neighbors and
K-Means algorithm and by the inclusion of an efficient algorithm for
computing the Ridge Generalized Cross Validation solution. Unlike the
preceding release, no new modules where added to this release.

Changelog
---------

  - Performance improvements for Gaussian Mixture Model sampling [Jan
    Schlüter].

  - Implementation of efficient leave-one-out cross-validated Ridge in
    :class:`linear_model.RidgeCV` [`Mathieu Blondel`_]

  - Better handling of collinearity and early stopping in
    :func:`linear_model.lars_path` [`Alexandre Gramfort`_ and `Fabian
    Pedregosa`_].

  - Fixes for liblinear ordering of labels and sign of coefficients
    [Dan Yamins, Paolo Losi, `Mathieu Blondel`_ and `Fabian Pedregosa`_].

  - Performance improvements for Nearest Neighbors algorithm in
    high-dimensional spaces [`Fabian Pedregosa`_].

  - Performance improvements for :class:`cluster.KMeans` [`Gael
    Varoquaux`_ and `James Bergstra`_].

  - Sanity checks for SVM-based classes [`Mathieu Blondel`_].

  - Refactoring of :class:`neighbors.NeighborsClassifier` and
    :func:`neighbors.kneighbors_graph`: added different algorithms for
    the k-Nearest Neighbor Search and implemented a more stable
    algorithm for finding barycenter weigths. Also added some
    developer documentation for this module, see
    :ref:`notes_neighbors` for more information [`Fabian Pedregosa`_].

  - Documentation improvements: Added :class:`pca.RandomizedPCA` and
    :class:`linear_model.LogisticRegression` to the class
    reference. Also added references of matrices used for clustering
    and other fixes [`Gael Varoquaux`_, `Fabian Pedregosa`_, `Mathieu
    Blondel`_, `Olivier Grisel`_, Virgile Fritsch , Emmanuelle
    Gouillart]

  - Binded decision_function in classes that make use of liblinear_,
    dense and sparse variants, like :class:`svm.LinearSVC` or
    :class:`linear_model.LogisticRegression` [`Fabian Pedregosa`_].

  - Performance and API improvements to
    :func:`metrics.euclidean_distances` and to
    :class:`pca.RandomizedPCA` [`James Bergstra`_].

  - Fix compilation issues under NetBSD [Kamel Ibn Hassen Derouiche]

  - Allow input sequences of different lengths in :class:`hmm.GaussianHMM`
    [`Ron Weiss`_].

  - Fix bug in affinity propagation caused by incorrect indexing [Xinfan Meng]


People
------

People that made this release possible preceeded by number of commits:

    - 85  `Fabian Pedregosa`_
    - 67  `Mathieu Blondel`_
    - 20  `Alexandre Gramfort`_
    - 19  `James Bergstra`_
    - 14  Dan Yamins
    - 13  `Olivier Grisel`_
    - 12  `Gael Varoquaux`_
    - 4  Edouard Duchesnay
    - 4  `Ron Weiss`_
    - 2  Satrajit Ghosh
    - 2  Vincent Dubourg
    - 1  Emmanuelle Gouillart
    - 1  Kamel Ibn Hassen Derouiche
    - 1  Paolo Losi
    - 1  VirgileFritsch
    - 1  `Yaroslav Halchenko`_
    - 1  Xinfan Meng


.. _changes_0_6:

0.6
===

scikits.learn 0.6 was released on december 2010. It is marked by the
inclusion of several new modules and a general renaming of old
ones. It is also marked by the inclusion of new example, including
applications to real-world datasets.

.. |banner1| image:: auto_examples/applications/images/plot_face_recognition_1.png
   :height: 150
   :target: auto_examples/applications/plot_face_recognition.html

.. |banner2| image:: auto_examples/applications/images/plot_species_distribution_modeling_1.png
   :height: 150
   :target: auto_examples/linear_model/plot_species_distribution.html

.. |banner3| image:: auto_examples/gaussian_process/images/plot_gp_regression_1.png
   :height: 150
   :target: auto_examples/gaussian_process/plot_gp_regression.html

.. |banner4| image:: auto_examples/linear_model/images/plot_sgd_iris_1.png
   :height: 150
   :target: auto_examples/linear_model/plot_lasso_lars.html


.. |center-div| raw:: html

    <div style="text-align: center; margin: 0px 0 -5px 0;">

.. |end-div| raw:: html

    </div>


|center-div| |banner1| |banner2| |banner3| |banner4| |end-div| 


Changelog
---------

  - New `stochastic gradient
    <http://scikit-learn.sourceforge.net/modules/sgd.html>`_ descent
    module by Peter Prettenhofer. The module comes with complete
    documentation and examples.

  - Improved svm module: memory consumption has been reduced by 50%,
    heuristic to automatically set class weights, possibility to
    assign weights to samples (see
    :ref:`example_svm_plot_weighted_samples.py` for an example).

  - New :ref:`gaussian_process` module by Vincent Dubourg. This module
    also has great documentation and some very neat examples. See
    :ref:`example_gaussian_process_plot_gp_regression.py` or
    :ref:`example_gaussian_process_plot_gp_probabilistic_classification_after_regression.py`
    for a taste of what can be done.

  - It is now possible to use liblinear’s Multi-class SVC (option
    multi_class in :class:`svm.LinearSVC`)

  - New features and performance improvements of text feature
    extraction. 

  - Improved sparse matrix support, both in main classes
    (:class:`grid_search.GridSearchCV`) as in modules
    scikits.learn.svm.sparse and scikits.learn.linear_model.sparse.

  - Lots of cool new examples and a new section that uses real-world
    datasets was created. These include:
    :ref:`example_applications_plot_face_recognition.py`,
    :ref:`example_applications_plot_species_distribution_modeling.py`,
    :ref:`example_applications_svm_gui.py`,
    :ref:`example_applications_wikipedia_principal_eigenvector.py` and
    others.

  - Faster :ref:`least_angle_regression` algorithm. It is now 2x
    faster than the R version on worst case and up to 10x times faster
    on some cases.

  - Faster coordinate descent algorithm. In particular, the full path
    version of lasso (:func:`linear_model.lasso_path`) is more than
    200x times faster than before.

  - It is now possible to get probability estimates from a
    :class:`linear_model.LogisticRegression` model.

  - module renaming: the glm module has been renamed to linear_model,
    the gmm module has been included into the more general mixture
    model and the sgd module has been included in linear_model.

  - Lots of bug fixes and documentation improvements.


People
------

People that made this release possible preceeded by number of commits:

   * 207  `Olivier Grisel <http://twitter.com/ogrisel>`_

   * 167 `Fabian Pedregosa <http://fseoane.net/blog/>`_

   * 97 `Peter Prettenhofer <http://sites.google.com/site/peterprettenhofer/>`_

   * 68 `Alexandre Gramfort <http://www-sop.inria.fr/members/Alexandre.Gramfort/>`_

   * 59  `Mathieu Blondel <http://www.mblondel.org/journal/>`_

   * 55  `Gael Varoquaux <http://gael-varoquaux.info/blog/>`_

   * 33  Vincent Dubourg

   * 21  `Ron Weiss <http://www.ee.columbia.edu/~ronw/>`_

   * 9  Bertrand Thirion

   * 3  `Alexandre Passos <http://atpassos.posterous.com>`_

   * 3  Anne-Laure Fouque

   * 2  Ronan Amicel

   * 1 `Christian Osendorfer <http://osdf.github.com/>`_



.. _changes_0_5:


0.5
===

Changelog
---------

New classes
~~~~~~~~~~~~

    - Support for sparse matrices in some classifiers of modules
      ``svm`` and ``linear_model`` (see :class:`svm.sparse.SVC`,
      :class:`svm.sparse.SVR`, :class:`svm.sparse.LinearSVC`,
      :class:`linear_model.sparse.Lasso`, :class:`linear_model.sparse.ElasticNet`)

    - New :class:`pipeline.Pipeline` object to compose different estimators.

    - Recursive Feature Elimination routines in module
      :ref:`feature_selection_doc`.

    - Addition of various classes capable of cross validation in the
      linear_model module (:class:`linear_model.LassoCV`, :class:`linear_model.ElasticNetCV`,
      etc.).

    - New, more efficient LARS algorithm implementation. The Lasso
      variant of the algorithm is also implemented. See
      :class:`linear_model.lars_path`, :class:`linear_model.LARS` and
      :class:`linear_model.LassoLARS`.

    - New Hidden Markov Models module (see classes
      :class:`hmm.GaussianHMM`, :class:`hmm.MultinomialHMM`,
      :class:`hmm.GMMHMM`)

    - New module feature_extraction (see :ref:`class reference
      <feature_extraction_ref>`)

    - New FastICA algorithm in module scikits.learn.fastica


Documentation
~~~~~~~~~~~~~

    - Improved documentation for many modules, now separating
      narrative documentation from the class reference. As an example,
      see `documentation for the SVM module
      <http://scikit-learn.sourceforge.net/modules/svm.html>`_ and the
      complete `class reference
      <http://scikit-learn.sourceforge.net/modules/classes.html>`_.

Fixes
~~~~~

    - API changes: adhere variable names to PEP-8, give more
      meaningful names.

    - Fixes for svm module to run on a shared memory context
      (multiprocessing).

    - It is again possible to generate latex (and thus PDF) from the
      sphinx docs.

Examples
~~~~~~~~

    - new examples using some of the mlcomp datasets:
      :ref:`example_mlcomp_sparse_document_classification.py`,
      :ref:`example_mlcomp_document_classification.py`

    - Many more examaples. `See here
      <http://scikit-learn.sourceforge.net/auto_examples/index.html>`_
      the full list of examples.


External dependencies
~~~~~~~~~~~~~~~~~~~~~

    - Joblib is now a dependencie of this package, although it is
      shipped with (scikits.learn.externals.joblib).

Removed modules
~~~~~~~~~~~~~~~

    - Module ann (Artificial Neural Networks) has been removed from
      the distribution. Users wanting this sort of algorithms should
      take a look into pybrain.

Misc
~~~~

    - New sphinx theme for the web page.


Authors
-------

The following is a list of authors for this release, preceeded by
number of commits:

     * 262  Fabian Pedregosa
     * 240  Gael Varoquaux
     * 149  Alexandre Gramfort
     * 116  Olivier Grisel
     *  40  Vincent Michel
     *  38  Ron Weiss
     *  23  Matthieu Perrot
     *  10  Bertrand Thirion
     *   7  Yaroslav Halchenko
     *   9  VirgileFritsch
     *   6  Edouard Duchesnay
     *   4  Mathieu Blondel
     *   1  Ariel Rokem
     *   1  Matthieu Brucher

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


.. _Alexandre Gramfort: http://www-sop.inria.fr/members/Alexandre.Gramfort/

.. _Fabian Pedregosa: http://fseoane.net/blog/

.. _Mathieu Blondel: http://www.mblondel.org/journal/

.. _James Bergstra: http://www-etud.iro.umontreal.ca/~bergstrj/

.. _liblinear: http://www.csie.ntu.edu.tw/~cjlin/liblinear/

.. _Yaroslav Halchenko: http://www.onerussian.com/
