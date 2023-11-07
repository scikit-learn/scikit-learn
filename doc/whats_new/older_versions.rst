.. include:: _contributors.rst

.. currentmodule:: sklearn

.. _changes_0_12.1:

Version 0.12.1
===============

**October 8, 2012**

The 0.12.1 release is a bug-fix release with no additional features, but is
instead a set of bug fixes

Changelog
----------

- Improved numerical stability in spectral embedding by `Gael
  Varoquaux`_

- Doctest under windows 64bit by `Gael Varoquaux`_

- Documentation fixes for elastic net by `Andreas Müller`_ and
  `Alexandre Gramfort`_

- Proper behavior with fortran-ordered NumPy arrays by `Gael Varoquaux`_

- Make GridSearchCV work with non-CSR sparse matrix by `Lars Buitinck`_

- Fix parallel computing in MDS by `Gael Varoquaux`_

- Fix Unicode support in count vectorizer by `Andreas Müller`_

- Fix MinCovDet breaking with X.shape = (3, 1) by :user:`Virgile Fritsch <VirgileFritsch>`

- Fix clone of SGD objects by `Peter Prettenhofer`_

- Stabilize GMM by :user:`Virgile Fritsch <VirgileFritsch>`

People
------

 *  14  `Peter Prettenhofer`_
 *  12  `Gael Varoquaux`_
 *  10  `Andreas Müller`_
 *   5  `Lars Buitinck`_
 *   3  :user:`Virgile Fritsch <VirgileFritsch>`
 *   1  `Alexandre Gramfort`_
 *   1  `Gilles Louppe`_
 *   1  `Mathieu Blondel`_

.. _changes_0_12:

Version 0.12
============

**September 4, 2012**

Changelog
---------

- Various speed improvements of the :ref:`decision trees <tree>` module, by
  `Gilles Louppe`_.

- :class:`~ensemble.GradientBoostingRegressor` and
  :class:`~ensemble.GradientBoostingClassifier` now support feature subsampling
  via the ``max_features`` argument, by `Peter Prettenhofer`_.

- Added Huber and Quantile loss functions to
  :class:`~ensemble.GradientBoostingRegressor`, by `Peter Prettenhofer`_.

- :ref:`Decision trees <tree>` and :ref:`forests of randomized trees <forest>`
  now support multi-output classification and regression problems, by
  `Gilles Louppe`_.

- Added :class:`~preprocessing.LabelEncoder`, a simple utility class to
  normalize labels or transform non-numerical labels, by `Mathieu Blondel`_.

- Added the epsilon-insensitive loss and the ability to make probabilistic
  predictions with the modified huber loss in :ref:`sgd`, by
  `Mathieu Blondel`_.

- Added :ref:`multidimensional_scaling`, by Nelle Varoquaux.

- SVMlight file format loader now detects compressed (gzip/bzip2) files and
  decompresses them on the fly, by `Lars Buitinck`_.

- SVMlight file format serializer now preserves double precision floating
  point values, by `Olivier Grisel`_.

- A common testing framework for all estimators was added, by `Andreas Müller`_.

- Understandable error messages for estimators that do not accept
  sparse input by `Gael Varoquaux`_

- Speedups in hierarchical clustering by `Gael Varoquaux`_. In
  particular building the tree now supports early stopping. This is
  useful when the number of clusters is not small compared to the
  number of samples.

- Add MultiTaskLasso and MultiTaskElasticNet for joint feature selection,
  by `Alexandre Gramfort`_.

- Added `metrics.auc_score` and
  :func:`metrics.average_precision_score` convenience functions by `Andreas
  Müller`_.

- Improved sparse matrix support in the :ref:`feature_selection`
  module by `Andreas Müller`_.

- New word boundaries-aware character n-gram analyzer for the
  :ref:`text_feature_extraction` module by :user:`@kernc <kernc>`.

- Fixed bug in spectral clustering that led to single point clusters
  by `Andreas Müller`_.

- In :class:`~feature_extraction.text.CountVectorizer`, added an option to
  ignore infrequent words, ``min_df`` by  `Andreas Müller`_.

- Add support for multiple targets in some linear models (ElasticNet, Lasso
  and OrthogonalMatchingPursuit) by `Vlad Niculae`_ and
  `Alexandre Gramfort`_.

- Fixes in `decomposition.ProbabilisticPCA` score function by Wei Li.

- Fixed feature importance computation in
  :ref:`gradient_boosting`.

API changes summary
-------------------

- The old ``scikits.learn`` package has disappeared; all code should import
  from ``sklearn`` instead, which was introduced in 0.9.

- In :func:`metrics.roc_curve`, the ``thresholds`` array is now returned
  with it's order reversed, in order to keep it consistent with the order
  of the returned ``fpr`` and ``tpr``.

- In `hmm` objects, like `hmm.GaussianHMM`,
  `hmm.MultinomialHMM`, etc., all parameters must be passed to the
  object when initialising it and not through ``fit``. Now ``fit`` will
  only accept the data as an input parameter.

- For all SVM classes, a faulty behavior of ``gamma`` was fixed. Previously,
  the default gamma value was only computed the first time ``fit`` was called
  and then stored. It is now recalculated on every call to ``fit``.

- All ``Base`` classes are now abstract meta classes so that they can not be
  instantiated.

- :func:`cluster.ward_tree` now also returns the parent array. This is
  necessary for early-stopping in which case the tree is not
  completely built.

- In :class:`~feature_extraction.text.CountVectorizer` the parameters
  ``min_n`` and ``max_n`` were joined to the parameter ``n_gram_range`` to
  enable grid-searching both at once.

- In :class:`~feature_extraction.text.CountVectorizer`, words that appear
  only in one document are now ignored by default. To reproduce
  the previous behavior, set ``min_df=1``.

- Fixed API inconsistency: :meth:`linear_model.SGDClassifier.predict_proba` now
  returns 2d array when fit on two classes.

- Fixed API inconsistency: :meth:`discriminant_analysis.QuadraticDiscriminantAnalysis.decision_function`
  and :meth:`discriminant_analysis.LinearDiscriminantAnalysis.decision_function` now return 1d arrays
  when fit on two classes.

- Grid of alphas used for fitting :class:`~linear_model.LassoCV` and
  :class:`~linear_model.ElasticNetCV` is now stored
  in the attribute ``alphas_`` rather than overriding the init parameter
  ``alphas``.

- Linear models when alpha is estimated by cross-validation store
  the estimated value in the ``alpha_`` attribute rather than just
  ``alpha`` or ``best_alpha``.

- :class:`~ensemble.GradientBoostingClassifier` now supports
  :meth:`~ensemble.GradientBoostingClassifier.staged_predict_proba`, and
  :meth:`~ensemble.GradientBoostingClassifier.staged_predict`.

- `svm.sparse.SVC` and other sparse SVM classes are now deprecated.
  The all classes in the :ref:`svm` module now automatically select the
  sparse or dense representation base on the input.

- All clustering algorithms now interpret the array ``X`` given to ``fit`` as
  input data, in particular :class:`~cluster.SpectralClustering` and
  :class:`~cluster.AffinityPropagation` which previously expected affinity matrices.

- For clustering algorithms that take the desired number of clusters as a parameter,
  this parameter is now called ``n_clusters``.


People
------
 * 267  `Andreas Müller`_
 *  94  `Gilles Louppe`_
 *  89  `Gael Varoquaux`_
 *  79  `Peter Prettenhofer`_
 *  60  `Mathieu Blondel`_
 *  57  `Alexandre Gramfort`_
 *  52  `Vlad Niculae`_
 *  45  `Lars Buitinck`_
 *  44  Nelle Varoquaux
 *  37  `Jaques Grobler`_
 *  30  Alexis Mignon
 *  30  Immanuel Bayer
 *  27  `Olivier Grisel`_
 *  16  Subhodeep Moitra
 *  13  Yannick Schwartz
 *  12  :user:`@kernc <kernc>`
 *  11  :user:`Virgile Fritsch <VirgileFritsch>`
 *   9  Daniel Duckworth
 *   9  `Fabian Pedregosa`_
 *   9  `Robert Layton`_
 *   8  John Benediktsson
 *   7  Marko Burjek
 *   5  `Nicolas Pinto`_
 *   4  Alexandre Abraham
 *   4  `Jake Vanderplas`_
 *   3  `Brian Holt`_
 *   3  `Edouard Duchesnay`_
 *   3  Florian Hoenig
 *   3  flyingimmidev
 *   2  Francois Savard
 *   2  Hannes Schulz
 *   2  Peter Welinder
 *   2  `Yaroslav Halchenko`_
 *   2  Wei Li
 *   1  Alex Companioni
 *   1  Brandyn A. White
 *   1  Bussonnier Matthias
 *   1  Charles-Pierre Astolfi
 *   1  Dan O'Huiginn
 *   1  David Cournapeau
 *   1  Keith Goodman
 *   1  Ludwig Schwardt
 *   1  Olivier Hervieu
 *   1  Sergio Medina
 *   1  Shiqiao Du
 *   1  Tim Sheerman-Chase
 *   1  buguen



.. _changes_0_11:

Version 0.11
============

**May 7, 2012**

Changelog
---------

Highlights
.............

- Gradient boosted regression trees (:ref:`gradient_boosting`)
  for classification and regression by `Peter Prettenhofer`_
  and `Scott White`_ .

- Simple dict-based feature loader with support for categorical variables
  (:class:`~feature_extraction.DictVectorizer`) by `Lars Buitinck`_.

- Added Matthews correlation coefficient (:func:`metrics.matthews_corrcoef`)
  and added macro and micro average options to
  :func:`~metrics.precision_score`, :func:`metrics.recall_score` and
  :func:`~metrics.f1_score` by `Satrajit Ghosh`_.

- :ref:`out_of_bag` of generalization error for :ref:`ensemble`
  by `Andreas Müller`_.

- Randomized sparse linear models for feature
  selection, by `Alexandre Gramfort`_ and `Gael Varoquaux`_

- :ref:`label_propagation` for semi-supervised learning, by Clay
  Woolam. **Note** the semi-supervised API is still work in progress,
  and may change.

- Added BIC/AIC model selection to classical :ref:`gmm` and unified
  the API with the remainder of scikit-learn, by `Bertrand Thirion`_

- Added `sklearn.cross_validation.StratifiedShuffleSplit`, which is
  a `sklearn.cross_validation.ShuffleSplit` with balanced splits,
  by Yannick Schwartz.

- :class:`~sklearn.neighbors.NearestCentroid` classifier added, along with a
  ``shrink_threshold`` parameter, which implements **shrunken centroid
  classification**, by `Robert Layton`_.

Other changes
..............

- Merged dense and sparse implementations of :ref:`sgd` module and
  exposed utility extension types for sequential
  datasets ``seq_dataset`` and weight vectors ``weight_vector``
  by `Peter Prettenhofer`_.

- Added ``partial_fit`` (support for online/minibatch learning) and
  warm_start to the :ref:`sgd` module by `Mathieu Blondel`_.

- Dense and sparse implementations of :ref:`svm` classes and
  :class:`~linear_model.LogisticRegression` merged by `Lars Buitinck`_.

- Regressors can now be used as base estimator in the :ref:`multiclass`
  module by `Mathieu Blondel`_.

- Added n_jobs option to :func:`metrics.pairwise_distances`
  and :func:`metrics.pairwise.pairwise_kernels` for parallel computation,
  by `Mathieu Blondel`_.

- :ref:`k_means` can now be run in parallel, using the ``n_jobs`` argument
  to either :ref:`k_means` or :class:`cluster.KMeans`, by `Robert Layton`_.

- Improved :ref:`cross_validation` and :ref:`grid_search` documentation
  and introduced the new `cross_validation.train_test_split`
  helper function by `Olivier Grisel`_

- :class:`~svm.SVC` members ``coef_`` and ``intercept_`` changed sign for
  consistency with ``decision_function``; for ``kernel==linear``,
  ``coef_`` was fixed in the one-vs-one case, by `Andreas Müller`_.

- Performance improvements to efficient leave-one-out cross-validated
  Ridge regression, esp. for the ``n_samples > n_features`` case, in
  :class:`~linear_model.RidgeCV`, by Reuben Fletcher-Costin.

- Refactoring and simplification of the :ref:`text_feature_extraction`
  API and fixed a bug that caused possible negative IDF,
  by `Olivier Grisel`_.

- Beam pruning option in `_BaseHMM` module has been removed since it
  is difficult to Cythonize. If you are interested in contributing a Cython
  version, you can use the python version in the git history as a reference.

- Classes in :ref:`neighbors` now support arbitrary Minkowski metric for
  nearest neighbors searches. The metric can be specified by argument ``p``.

API changes summary
-------------------

- `covariance.EllipticEnvelop` is now deprecated.
  Please use :class:`~covariance.EllipticEnvelope` instead.

- ``NeighborsClassifier`` and ``NeighborsRegressor`` are gone in the module
  :ref:`neighbors`. Use the classes :class:`~neighbors.KNeighborsClassifier`,
  :class:`~neighbors.RadiusNeighborsClassifier`, :class:`~neighbors.KNeighborsRegressor`
  and/or :class:`~neighbors.RadiusNeighborsRegressor` instead.

- Sparse classes in the :ref:`sgd` module are now deprecated.

- In `mixture.GMM`, `mixture.DPGMM` and `mixture.VBGMM`,
  parameters must be passed to an object when initialising it and not through
  ``fit``. Now ``fit`` will only accept the data as an input parameter.

- methods ``rvs`` and ``decode`` in `GMM` module are now deprecated.
  ``sample`` and ``score`` or ``predict`` should be used instead.

- attribute ``_scores`` and ``_pvalues`` in univariate feature selection
  objects are now deprecated.
  ``scores_`` or ``pvalues_`` should be used instead.

- In :class:`~linear_model.LogisticRegression`, :class:`~svm.LinearSVC`,
  :class:`~svm.SVC` and :class:`~svm.NuSVC`, the ``class_weight`` parameter is
  now an initialization parameter, not a parameter to fit. This makes grid
  searches over this parameter possible.

- LFW ``data`` is now always shape ``(n_samples, n_features)`` to be
  consistent with the Olivetti faces dataset. Use ``images`` and
  ``pairs`` attribute to access the natural images shapes instead.

- In :class:`~svm.LinearSVC`, the meaning of the ``multi_class`` parameter
  changed.  Options now are ``'ovr'`` and ``'crammer_singer'``, with
  ``'ovr'`` being the default.  This does not change the default behavior
  but hopefully is less confusing.

- Class `feature_selection.text.Vectorizer` is deprecated and
  replaced by `feature_selection.text.TfidfVectorizer`.

- The preprocessor / analyzer nested structure for text feature
  extraction has been removed. All those features are
  now directly passed as flat constructor arguments
  to `feature_selection.text.TfidfVectorizer` and
  `feature_selection.text.CountVectorizer`, in particular the
  following parameters are now used:

- ``analyzer`` can be ``'word'`` or ``'char'`` to switch the default
  analysis scheme, or use a specific python callable (as previously).

- ``tokenizer`` and ``preprocessor`` have been introduced to make it
  still possible to customize those steps with the new API.

- ``input`` explicitly control how to interpret the sequence passed to
  ``fit`` and ``predict``: filenames, file objects or direct (byte or
  Unicode) strings.

- charset decoding is explicit and strict by default.

- the ``vocabulary``, fitted or not is now stored in the
  ``vocabulary_`` attribute to be consistent with the project
  conventions.

- Class `feature_selection.text.TfidfVectorizer` now derives directly
  from `feature_selection.text.CountVectorizer` to make grid
  search trivial.

- methods ``rvs`` in `_BaseHMM` module are now deprecated.
  ``sample`` should be used instead.

- Beam pruning option in `_BaseHMM` module is removed since it is
  difficult to be Cythonized. If you are interested, you can look in the
  history codes by git.

- The SVMlight format loader now supports files with both zero-based and
  one-based column indices, since both occur "in the wild".

- Arguments in class :class:`~model_selection.ShuffleSplit` are now consistent with
  :class:`~model_selection.StratifiedShuffleSplit`. Arguments ``test_fraction`` and
  ``train_fraction`` are deprecated and renamed to ``test_size`` and
  ``train_size`` and can accept both ``float`` and ``int``.

- Arguments in class `Bootstrap` are now consistent with
  :class:`~model_selection.StratifiedShuffleSplit`. Arguments ``n_test`` and
  ``n_train`` are deprecated and renamed to ``test_size`` and
  ``train_size`` and can accept both ``float`` and ``int``.

- Argument ``p`` added to classes in :ref:`neighbors` to specify an
  arbitrary Minkowski metric for nearest neighbors searches.


People
------
   * 282  `Andreas Müller`_
   * 239  `Peter Prettenhofer`_
   * 198  `Gael Varoquaux`_
   * 129  `Olivier Grisel`_
   * 114  `Mathieu Blondel`_
   * 103  Clay Woolam
   *  96  `Lars Buitinck`_
   *  88  `Jaques Grobler`_
   *  82  `Alexandre Gramfort`_
   *  50  `Bertrand Thirion`_
   *  42  `Robert Layton`_
   *  28  flyingimmidev
   *  26  `Jake Vanderplas`_
   *  26  Shiqiao Du
   *  21  `Satrajit Ghosh`_
   *  17  `David Marek`_
   *  17  `Gilles Louppe`_
   *  14  `Vlad Niculae`_
   *  11  Yannick Schwartz
   *  10  `Fabian Pedregosa`_
   *   9  fcostin
   *   7  Nick Wilson
   *   5  Adrien Gaidon
   *   5  `Nicolas Pinto`_
   *   4  `David Warde-Farley`_
   *   5  Nelle Varoquaux
   *   5  Emmanuelle Gouillart
   *   3  Joonas Sillanpää
   *   3  Paolo Losi
   *   2  Charles McCarthy
   *   2  Roy Hyunjin Han
   *   2  Scott White
   *   2  ibayer
   *   1  Brandyn White
   *   1  Carlos Scheidegger
   *   1  Claire Revillet
   *   1  Conrad Lee
   *   1  `Edouard Duchesnay`_
   *   1  Jan Hendrik Metzen
   *   1  Meng Xinfan
   *   1  `Rob Zinkov`_
   *   1  Shiqiao
   *   1  Udi Weinsberg
   *   1  Virgile Fritsch
   *   1  Xinfan Meng
   *   1  Yaroslav Halchenko
   *   1  jansoe
   *   1  Leon Palafox


.. _changes_0_10:

Version 0.10
============

**January 11, 2012**

Changelog
---------

- Python 2.5 compatibility was dropped; the minimum Python version needed
  to use scikit-learn is now 2.6.

- :ref:`sparse_inverse_covariance` estimation using the graph Lasso, with
  associated cross-validated estimator, by `Gael Varoquaux`_

- New :ref:`Tree <tree>` module by `Brian Holt`_, `Peter Prettenhofer`_,
  `Satrajit Ghosh`_ and `Gilles Louppe`_. The module comes with complete
  documentation and examples.

- Fixed a bug in the RFE module by `Gilles Louppe`_ (issue #378).

- Fixed a memory leak in :ref:`svm` module by `Brian Holt`_ (issue #367).

- Faster tests by `Fabian Pedregosa`_ and others.

- Silhouette Coefficient cluster analysis evaluation metric added as
  :func:`~sklearn.metrics.silhouette_score` by Robert Layton.

- Fixed a bug in :ref:`k_means` in the handling of the ``n_init`` parameter:
  the clustering algorithm used to be run ``n_init`` times but the last
  solution was retained instead of the best solution by `Olivier Grisel`_.

- Minor refactoring in :ref:`sgd` module; consolidated dense and sparse
  predict methods; Enhanced test time performance by converting model
  parameters to fortran-style arrays after fitting (only multi-class).

- Adjusted Mutual Information metric added as
  :func:`~sklearn.metrics.adjusted_mutual_info_score` by Robert Layton.

- Models like SVC/SVR/LinearSVC/LogisticRegression from libsvm/liblinear
  now support scaling of C regularization parameter by the number of
  samples by `Alexandre Gramfort`_.

- New :ref:`Ensemble Methods <ensemble>` module by `Gilles Louppe`_ and
  `Brian Holt`_. The module comes with the random forest algorithm and the
  extra-trees method, along with documentation and examples.

- :ref:`outlier_detection`: outlier and novelty detection, by
  :user:`Virgile Fritsch <VirgileFritsch>`.

- :ref:`kernel_approximation`: a transform implementing kernel
  approximation for fast SGD on non-linear kernels by
  `Andreas Müller`_.

- Fixed a bug due to atom swapping in :ref:`OMP` by `Vlad Niculae`_.

- :ref:`SparseCoder` by `Vlad Niculae`_.

- :ref:`mini_batch_kmeans` performance improvements by `Olivier Grisel`_.

- :ref:`k_means` support for sparse matrices by `Mathieu Blondel`_.

- Improved documentation for developers and for the :mod:`sklearn.utils`
  module, by `Jake Vanderplas`_.

- Vectorized 20newsgroups dataset loader
  (:func:`~sklearn.datasets.fetch_20newsgroups_vectorized`) by
  `Mathieu Blondel`_.

- :ref:`multiclass` by `Lars Buitinck`_.

- Utilities for fast computation of mean and variance for sparse matrices
  by `Mathieu Blondel`_.

- Make :func:`~sklearn.preprocessing.scale` and
  `sklearn.preprocessing.Scaler` work on sparse matrices by
  `Olivier Grisel`_

- Feature importances using decision trees and/or forest of trees,
  by `Gilles Louppe`_.

- Parallel implementation of forests of randomized trees by
  `Gilles Louppe`_.

- `sklearn.cross_validation.ShuffleSplit` can subsample the train
  sets as well as the test sets by `Olivier Grisel`_.

- Errors in the build of the documentation fixed by `Andreas Müller`_.


API changes summary
-------------------

Here are the code migration instructions when upgrading from scikit-learn
version 0.9:

- Some estimators that may overwrite their inputs to save memory previously
  had ``overwrite_`` parameters; these have been replaced with ``copy_``
  parameters with exactly the opposite meaning.

  This particularly affects some of the estimators in :mod:`~sklearn.linear_model`.
  The default behavior is still to copy everything passed in.

- The SVMlight dataset loader :func:`~sklearn.datasets.load_svmlight_file` no
  longer supports loading two files at once; use ``load_svmlight_files``
  instead. Also, the (unused) ``buffer_mb`` parameter is gone.

- Sparse estimators in the :ref:`sgd` module use dense parameter vector
  ``coef_`` instead of ``sparse_coef_``. This significantly improves
  test time performance.

- The :ref:`covariance` module now has a robust estimator of
  covariance, the Minimum Covariance Determinant estimator.

- Cluster evaluation metrics in :mod:`~sklearn.metrics.cluster` have been refactored
  but the changes are backwards compatible. They have been moved to the
  `metrics.cluster.supervised`, along with
  `metrics.cluster.unsupervised` which contains the Silhouette
  Coefficient.

- The ``permutation_test_score`` function now behaves the same way as
  ``cross_val_score`` (i.e. uses the mean score across the folds.)

- Cross Validation generators now use integer indices (``indices=True``)
  by default instead of boolean masks. This make it more intuitive to
  use with sparse matrix data.

- The functions used for sparse coding, ``sparse_encode`` and
  ``sparse_encode_parallel`` have been combined into
  :func:`~sklearn.decomposition.sparse_encode`, and the shapes of the arrays
  have been transposed for consistency with the matrix factorization setting,
  as opposed to the regression setting.

- Fixed an off-by-one error in the SVMlight/LibSVM file format handling;
  files generated using :func:`~sklearn.datasets.dump_svmlight_file` should be
  re-generated. (They should continue to work, but accidentally had one
  extra column of zeros prepended.)

- ``BaseDictionaryLearning`` class replaced by ``SparseCodingMixin``.

- `sklearn.utils.extmath.fast_svd` has been renamed
  :func:`~sklearn.utils.extmath.randomized_svd` and the default
  oversampling is now fixed to 10 additional random vectors instead
  of doubling the number of components to extract. The new behavior
  follows the reference paper.


People
------

The following people contributed to scikit-learn since last release:

   * 246  `Andreas Müller`_
   * 242  `Olivier Grisel`_
   * 220  `Gilles Louppe`_
   * 183  `Brian Holt`_
   * 166  `Gael Varoquaux`_
   * 144  `Lars Buitinck`_
   *  73  `Vlad Niculae`_
   *  65  `Peter Prettenhofer`_
   *  64  `Fabian Pedregosa`_
   *  60  Robert Layton
   *  55  `Mathieu Blondel`_
   *  52  `Jake Vanderplas`_
   *  44  Noel Dawe
   *  38  `Alexandre Gramfort`_
   *  24  :user:`Virgile Fritsch <VirgileFritsch>`
   *  23  `Satrajit Ghosh`_
   *   3  Jan Hendrik Metzen
   *   3  Kenneth C. Arnold
   *   3  Shiqiao Du
   *   3  Tim Sheerman-Chase
   *   3  `Yaroslav Halchenko`_
   *   2  Bala Subrahmanyam Varanasi
   *   2  DraXus
   *   2  Michael Eickenberg
   *   1  Bogdan Trach
   *   1  Félix-Antoine Fortin
   *   1  Juan Manuel Caicedo Carvajal
   *   1  Nelle Varoquaux
   *   1  `Nicolas Pinto`_
   *   1  Tiziano Zito
   *   1  Xinfan Meng



.. _changes_0_9:

Version 0.9
===========

**September 21, 2011**

scikit-learn 0.9 was released on September 2011, three months after the 0.8
release and includes the new modules :ref:`manifold`, :ref:`dirichlet_process`
as well as several new algorithms and documentation improvements.

This release also includes the dictionary-learning work developed by
`Vlad Niculae`_ as part of the `Google Summer of Code
<https://developers.google.com/open-source/gsoc>`_ program.



.. |banner1| image:: ../auto_examples/manifold/images/thumb/sphx_glr_plot_compare_methods_thumb.png
   :target: ../auto_examples/manifold/plot_compare_methods.html

.. |banner2| image:: ../auto_examples/linear_model/images/thumb/sphx_glr_plot_omp_thumb.png
   :target: ../auto_examples/linear_model/plot_omp.html

.. |banner3| image:: ../auto_examples/decomposition/images/thumb/sphx_glr_plot_kernel_pca_thumb.png
   :target: ../auto_examples/decomposition/plot_kernel_pca.html

.. |center-div| raw:: html

    <div style="text-align: center; margin: 0px 0 -5px 0;">

.. |end-div| raw:: html

    </div>


|center-div| |banner2| |banner1| |banner3| |end-div|

Changelog
---------

- New :ref:`manifold` module by `Jake Vanderplas`_ and
  `Fabian Pedregosa`_.

- New :ref:`Dirichlet Process <dirichlet_process>` Gaussian Mixture
  Model by `Alexandre Passos`_

- :ref:`neighbors` module refactoring by `Jake Vanderplas`_ :
  general refactoring, support for sparse matrices in input, speed and
  documentation improvements. See the next section for a full list of API
  changes.

- Improvements on the :ref:`feature_selection` module by
  `Gilles Louppe`_ : refactoring of the RFE classes, documentation
  rewrite, increased efficiency and minor API changes.

- :ref:`SparsePCA` by `Vlad Niculae`_, `Gael Varoquaux`_ and
  `Alexandre Gramfort`_

- Printing an estimator now behaves independently of architectures
  and Python version thanks to :user:`Jean Kossaifi <JeanKossaifi>`.

- :ref:`Loader for libsvm/svmlight format <libsvm_loader>` by
  `Mathieu Blondel`_ and `Lars Buitinck`_

- Documentation improvements: thumbnails in
  example gallery by `Fabian Pedregosa`_.

- Important bugfixes in :ref:`svm` module (segfaults, bad
  performance) by `Fabian Pedregosa`_.

- Added :ref:`multinomial_naive_bayes` and :ref:`bernoulli_naive_bayes`
  by `Lars Buitinck`_

- Text feature extraction optimizations by Lars Buitinck

- Chi-Square feature selection
  (:func:`feature_selection.chi2`) by `Lars Buitinck`_.

- :ref:`sample_generators` module refactoring by `Gilles Louppe`_

- :ref:`multiclass` by `Mathieu Blondel`_

- Ball tree rewrite by `Jake Vanderplas`_

- Implementation of :ref:`dbscan` algorithm by Robert Layton

- Kmeans predict and transform by Robert Layton

- Preprocessing module refactoring by `Olivier Grisel`_

- Faster mean shift by Conrad Lee

- New ``Bootstrap``, :ref:`ShuffleSplit` and various other
  improvements in cross validation schemes by `Olivier Grisel`_ and
  `Gael Varoquaux`_

- Adjusted Rand index and V-Measure clustering evaluation metrics by `Olivier Grisel`_

- Added :class:`Orthogonal Matching Pursuit <linear_model.OrthogonalMatchingPursuit>` by `Vlad Niculae`_

- Added 2D-patch extractor utilities in the :ref:`feature_extraction` module by `Vlad Niculae`_

- Implementation of :class:`~linear_model.LassoLarsCV`
  (cross-validated Lasso solver using the Lars algorithm) and
  :class:`~linear_model.LassoLarsIC` (BIC/AIC model
  selection in Lars) by `Gael Varoquaux`_
  and `Alexandre Gramfort`_

- Scalability improvements to :func:`metrics.roc_curve` by Olivier Hervieu

- Distance helper functions :func:`metrics.pairwise_distances`
  and :func:`metrics.pairwise.pairwise_kernels` by Robert Layton

- :class:`Mini-Batch K-Means <cluster.MiniBatchKMeans>` by Nelle Varoquaux and Peter Prettenhofer.

- mldata utilities by Pietro Berkes.

- :ref:`olivetti_faces_dataset` by `David Warde-Farley`_.


API changes summary
-------------------

Here are the code migration instructions when upgrading from scikit-learn
version 0.8:

- The ``scikits.learn`` package was renamed ``sklearn``. There is
  still a ``scikits.learn`` package alias for backward compatibility.

  Third-party projects with a dependency on scikit-learn 0.9+ should
  upgrade their codebase. For instance, under Linux / MacOSX just run
  (make a backup first!)::

      find -name "*.py" | xargs sed -i 's/\bscikits.learn\b/sklearn/g'

- Estimators no longer accept model parameters as ``fit`` arguments:
  instead all parameters must be only be passed as constructor
  arguments or using the now public ``set_params`` method inherited
  from :class:`~base.BaseEstimator`.

  Some estimators can still accept keyword arguments on the ``fit``
  but this is restricted to data-dependent values (e.g. a Gram matrix
  or an affinity matrix that are precomputed from the ``X`` data matrix.

- The ``cross_val`` package has been renamed to ``cross_validation``
  although there is also a ``cross_val`` package alias in place for
  backward compatibility.

  Third-party projects with a dependency on scikit-learn 0.9+ should
  upgrade their codebase. For instance, under Linux / MacOSX just run
  (make a backup first!)::

      find -name "*.py" | xargs sed -i 's/\bcross_val\b/cross_validation/g'

- The ``score_func`` argument of the
  ``sklearn.cross_validation.cross_val_score`` function is now expected
  to accept ``y_test`` and ``y_predicted`` as only arguments for
  classification and regression tasks or ``X_test`` for unsupervised
  estimators.

- ``gamma`` parameter for support vector machine algorithms is set
  to ``1 / n_features`` by default, instead of ``1 / n_samples``.

- The ``sklearn.hmm`` has been marked as orphaned: it will be removed
  from scikit-learn in version 0.11 unless someone steps up to
  contribute documentation, examples and fix lurking numerical
  stability issues.

- ``sklearn.neighbors`` has been made into a submodule.  The two previously
  available estimators, ``NeighborsClassifier`` and ``NeighborsRegressor``
  have been marked as deprecated.  Their functionality has been divided
  among five new classes: ``NearestNeighbors`` for unsupervised neighbors
  searches, ``KNeighborsClassifier`` & ``RadiusNeighborsClassifier``
  for supervised classification problems, and ``KNeighborsRegressor``
  & ``RadiusNeighborsRegressor`` for supervised regression problems.

- ``sklearn.ball_tree.BallTree`` has been moved to
  ``sklearn.neighbors.BallTree``.  Using the former will generate a warning.

- ``sklearn.linear_model.LARS()`` and related classes (LassoLARS,
  LassoLARSCV, etc.) have been renamed to
  ``sklearn.linear_model.Lars()``.

- All distance metrics and kernels in ``sklearn.metrics.pairwise`` now have a Y
  parameter, which by default is None. If not given, the result is the distance
  (or kernel similarity) between each sample in Y. If given, the result is the
  pairwise distance (or kernel similarity) between samples in X to Y.

- ``sklearn.metrics.pairwise.l1_distance`` is now called ``manhattan_distance``,
  and by default returns the pairwise distance. For the component wise distance,
  set the parameter ``sum_over_features`` to ``False``.

Backward compatibility package aliases and other deprecated classes and
functions will be removed in version 0.11.


People
------

38 people contributed to this release.

- 387  `Vlad Niculae`_
- 320  `Olivier Grisel`_
- 192  `Lars Buitinck`_
- 179  `Gael Varoquaux`_
- 168  `Fabian Pedregosa`_ (`INRIA`_, `Parietal Team`_)
- 127  `Jake Vanderplas`_
- 120  `Mathieu Blondel`_
- 85  `Alexandre Passos`_
- 67  `Alexandre Gramfort`_
- 57  `Peter Prettenhofer`_
- 56  `Gilles Louppe`_
- 42  Robert Layton
- 38  Nelle Varoquaux
- 32  :user:`Jean Kossaifi <JeanKossaifi>`
- 30  Conrad Lee
- 22  Pietro Berkes
- 18  andy
- 17  David Warde-Farley
- 12  Brian Holt
- 11  Robert
- 8  Amit Aides
- 8  :user:`Virgile Fritsch <VirgileFritsch>`
- 7  `Yaroslav Halchenko`_
- 6  Salvatore Masecchia
- 5  Paolo Losi
- 4  Vincent Schut
- 3  Alexis Metaireau
- 3  Bryan Silverthorn
- 3  `Andreas Müller`_
- 2  Minwoo Jake Lee
- 1  Emmanuelle Gouillart
- 1  Keith Goodman
- 1  Lucas Wiman
- 1  `Nicolas Pinto`_
- 1  Thouis (Ray) Jones
- 1  Tim Sheerman-Chase


.. _changes_0_8:

Version 0.8
===========

**May 11, 2011**

scikit-learn 0.8 was released on May 2011, one month after the first
"international" `scikit-learn coding sprint
<https://github.com/scikit-learn/scikit-learn/wiki/Upcoming-events>`_ and is
marked by the inclusion of important modules: :ref:`hierarchical_clustering`,
:ref:`cross_decomposition`, :ref:`NMF`, initial support for Python 3 and by important
enhancements and bug fixes.


Changelog
---------

Several new modules where introduced during this release:

- New :ref:`hierarchical_clustering` module by Vincent Michel,
  `Bertrand Thirion`_, `Alexandre Gramfort`_ and `Gael Varoquaux`_.

- :ref:`kernel_pca` implementation by `Mathieu Blondel`_

- :ref:`labeled_faces_in_the_wild_dataset` by `Olivier Grisel`_.

- New :ref:`cross_decomposition` module by `Edouard Duchesnay`_.

- :ref:`NMF` module `Vlad Niculae`_

- Implementation of the :ref:`oracle_approximating_shrinkage` algorithm by
  :user:`Virgile Fritsch <VirgileFritsch>` in the :ref:`covariance` module.


Some other modules benefited from significant improvements or cleanups.


- Initial support for Python 3: builds and imports cleanly,
  some modules are usable while others have failing tests by `Fabian Pedregosa`_.

- :class:`~decomposition.PCA` is now usable from the Pipeline object by `Olivier Grisel`_.

- Guide :ref:`performance-howto` by `Olivier Grisel`_.

- Fixes for memory leaks in libsvm bindings, 64-bit safer BallTree by Lars Buitinck.

- bug and style fixing in :ref:`k_means` algorithm by Jan Schlüter.

- Add attribute converged to Gaussian Mixture Models by Vincent Schut.

- Implemented ``transform``, ``predict_log_proba`` in
  :class:`~discriminant_analysis.LinearDiscriminantAnalysis` By `Mathieu Blondel`_.

- Refactoring in the :ref:`svm` module and bug fixes by `Fabian Pedregosa`_,
  `Gael Varoquaux`_ and Amit Aides.

- Refactored SGD module (removed code duplication, better variable naming),
  added interface for sample weight by `Peter Prettenhofer`_.

- Wrapped BallTree with Cython by Thouis (Ray) Jones.

- Added function :func:`svm.l1_min_c` by Paolo Losi.

- Typos, doc style, etc. by `Yaroslav Halchenko`_, `Gael Varoquaux`_,
  `Olivier Grisel`_, Yann Malet, `Nicolas Pinto`_, Lars Buitinck and
  `Fabian Pedregosa`_.


People
-------

People that made this release possible preceded by number of commits:


- 159  `Olivier Grisel`_
- 96  `Gael Varoquaux`_
- 96  `Vlad Niculae`_
- 94  `Fabian Pedregosa`_
- 36  `Alexandre Gramfort`_
- 32  Paolo Losi
- 31  `Edouard Duchesnay`_
- 30  `Mathieu Blondel`_
- 25  `Peter Prettenhofer`_
- 22  `Nicolas Pinto`_
- 11  :user:`Virgile Fritsch <VirgileFritsch>`
   -  7  Lars Buitinck
   -  6  Vincent Michel
   -  5  `Bertrand Thirion`_
   -  4  Thouis (Ray) Jones
   -  4  Vincent Schut
   -  3  Jan Schlüter
   -  2  Julien Miotte
   -  2  `Matthieu Perrot`_
   -  2  Yann Malet
   -  2  `Yaroslav Halchenko`_
   -  1  Amit Aides
   -  1  `Andreas Müller`_
   -  1  Feth Arezki
   -  1  Meng Xinfan


.. _changes_0_7:

Version 0.7
===========

**March 2, 2011**

scikit-learn 0.7 was released in March 2011, roughly three months
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
  :class:`~linear_model.RidgeCV` [`Mathieu Blondel`_]

- Better handling of collinearity and early stopping in
  :func:`linear_model.lars_path` [`Alexandre Gramfort`_ and `Fabian
  Pedregosa`_].

- Fixes for liblinear ordering of labels and sign of coefficients
  [Dan Yamins, Paolo Losi, `Mathieu Blondel`_ and `Fabian Pedregosa`_].

- Performance improvements for Nearest Neighbors algorithm in
  high-dimensional spaces [`Fabian Pedregosa`_].

- Performance improvements for :class:`~cluster.KMeans` [`Gael
  Varoquaux`_ and `James Bergstra`_].

- Sanity checks for SVM-based classes [`Mathieu Blondel`_].

- Refactoring of `neighbors.NeighborsClassifier` and
  :func:`neighbors.kneighbors_graph`: added different algorithms for
  the k-Nearest Neighbor Search and implemented a more stable
  algorithm for finding barycenter weights. Also added some
  developer documentation for this module, see
  `notes_neighbors
  <https://github.com/scikit-learn/scikit-learn/wiki/Neighbors-working-notes>`_ for more information [`Fabian Pedregosa`_].

- Documentation improvements: Added `pca.RandomizedPCA` and
  :class:`~linear_model.LogisticRegression` to the class
  reference. Also added references of matrices used for clustering
  and other fixes [`Gael Varoquaux`_, `Fabian Pedregosa`_, `Mathieu
  Blondel`_, `Olivier Grisel`_, Virgile Fritsch , Emmanuelle
  Gouillart]

- Binded decision_function in classes that make use of liblinear_,
  dense and sparse variants, like :class:`~svm.LinearSVC` or
  :class:`~linear_model.LogisticRegression` [`Fabian Pedregosa`_].

- Performance and API improvements to
  :func:`metrics.pairwise.euclidean_distances` and to
  `pca.RandomizedPCA` [`James Bergstra`_].

- Fix compilation issues under NetBSD [Kamel Ibn Hassen Derouiche]

- Allow input sequences of different lengths in `hmm.GaussianHMM`
  [`Ron Weiss`_].

- Fix bug in affinity propagation caused by incorrect indexing [Xinfan Meng]


People
------

People that made this release possible preceded by number of commits:

- 85  `Fabian Pedregosa`_
- 67  `Mathieu Blondel`_
- 20  `Alexandre Gramfort`_
- 19  `James Bergstra`_
- 14  Dan Yamins
- 13  `Olivier Grisel`_
- 12  `Gael Varoquaux`_
- 4  `Edouard Duchesnay`_
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

Version 0.6
===========

**December 21, 2010**

scikit-learn 0.6 was released on December 2010. It is marked by the
inclusion of several new modules and a general renaming of old
ones. It is also marked by the inclusion of new example, including
applications to real-world datasets.


Changelog
---------

- New `stochastic gradient
  <http://scikit-learn.org/stable/modules/sgd.html>`_ descent
  module by Peter Prettenhofer. The module comes with complete
  documentation and examples.

- Improved svm module: memory consumption has been reduced by 50%,
  heuristic to automatically set class weights, possibility to
  assign weights to samples (see
  :ref:`sphx_glr_auto_examples_svm_plot_weighted_samples.py` for an example).

- New :ref:`gaussian_process` module by Vincent Dubourg. This module
  also has great documentation and some very neat examples. See
  example_gaussian_process_plot_gp_regression.py or
  example_gaussian_process_plot_gp_probabilistic_classification_after_regression.py
  for a taste of what can be done.

- It is now possible to use liblinear's Multi-class SVC (option
  multi_class in :class:`~svm.LinearSVC`)

- New features and performance improvements of text feature
  extraction.

- Improved sparse matrix support, both in main classes
  (:class:`~model_selection.GridSearchCV`) as in modules
  sklearn.svm.sparse and sklearn.linear_model.sparse.

- Lots of cool new examples and a new section that uses real-world
  datasets was created. These include:
  :ref:`sphx_glr_auto_examples_applications_plot_face_recognition.py`,
  :ref:`sphx_glr_auto_examples_applications_plot_species_distribution_modeling.py`,
  :ref:`sphx_glr_auto_examples_applications_svm_gui.py`,
  :ref:`sphx_glr_auto_examples_applications_wikipedia_principal_eigenvector.py` and
  others.

- Faster :ref:`least_angle_regression` algorithm. It is now 2x
  faster than the R version on worst case and up to 10x times faster
  on some cases.

- Faster coordinate descent algorithm. In particular, the full path
  version of lasso (:func:`linear_model.lasso_path`) is more than
  200x times faster than before.

- It is now possible to get probability estimates from a
  :class:`~linear_model.LogisticRegression` model.

- module renaming: the glm module has been renamed to linear_model,
  the gmm module has been included into the more general mixture
  model and the sgd module has been included in linear_model.

- Lots of bug fixes and documentation improvements.


People
------

People that made this release possible preceded by number of commits:

   * 207  `Olivier Grisel`_

   * 167 `Fabian Pedregosa`_

   * 97 `Peter Prettenhofer`_

   * 68 `Alexandre Gramfort`_

   * 59  `Mathieu Blondel`_

   * 55  `Gael Varoquaux`_

   * 33  Vincent Dubourg

   * 21  `Ron Weiss`_

   * 9  Bertrand Thirion

   * 3  `Alexandre Passos`_

   * 3  Anne-Laure Fouque

   * 2  Ronan Amicel

   * 1 `Christian Osendorfer`_



.. _changes_0_5:


Version 0.5
===========

**October 11, 2010**

Changelog
---------

New classes
-----------

- Support for sparse matrices in some classifiers of modules
  ``svm`` and ``linear_model`` (see `svm.sparse.SVC`,
  `svm.sparse.SVR`, `svm.sparse.LinearSVC`,
  `linear_model.sparse.Lasso`, `linear_model.sparse.ElasticNet`)

- New :class:`~pipeline.Pipeline` object to compose different estimators.

- Recursive Feature Elimination routines in module
  :ref:`feature_selection`.

- Addition of various classes capable of cross validation in the
  linear_model module (:class:`~linear_model.LassoCV`, :class:`~linear_model.ElasticNetCV`,
  etc.).

- New, more efficient LARS algorithm implementation. The Lasso
  variant of the algorithm is also implemented. See
  :class:`~linear_model.lars_path`, :class:`~linear_model.Lars` and
  :class:`~linear_model.LassoLars`.

- New Hidden Markov Models module (see classes
  `hmm.GaussianHMM`, `hmm.MultinomialHMM`, `hmm.GMMHMM`)

- New module feature_extraction (see :ref:`class reference
  <feature_extraction_ref>`)

- New FastICA algorithm in module sklearn.fastica


Documentation
-------------

- Improved documentation for many modules, now separating
  narrative documentation from the class reference. As an example,
  see `documentation for the SVM module
  <http://scikit-learn.org/stable/modules/svm.html>`_ and the
  complete `class reference
  <http://scikit-learn.org/stable/modules/classes.html>`_.

Fixes
-----

- API changes: adhere variable names to PEP-8, give more
  meaningful names.

- Fixes for svm module to run on a shared memory context
  (multiprocessing).

- It is again possible to generate latex (and thus PDF) from the
  sphinx docs.

Examples
--------

- new examples using some of the mlcomp datasets:
  ``sphx_glr_auto_examples_mlcomp_sparse_document_classification.py`` (since removed) and
  :ref:`sphx_glr_auto_examples_text_plot_document_classification_20newsgroups.py`

- Many more examples. `See here
  <http://scikit-learn.org/stable/auto_examples/index.html>`_
  the full list of examples.


External dependencies
---------------------

- Joblib is now a dependency of this package, although it is
  shipped with (sklearn.externals.joblib).

Removed modules
---------------

- Module ann (Artificial Neural Networks) has been removed from
  the distribution. Users wanting this sort of algorithms should
  take a look into pybrain.

Misc
----

- New sphinx theme for the web page.


Authors
-------

The following is a list of authors for this release, preceded by
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

Version 0.4
===========

**August 26, 2010**

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

- Migration to GIT as version control system.

- Removal of obsolete attrselect module.

- Rename of private compiled extensions (added underscore).

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


Earlier versions
================

Earlier versions included contributions by Fred Mailhot, David Cooke,
David Huard, Dave Morrill, Ed Schofield, Travis Oliphant, Pearu Peterson.
