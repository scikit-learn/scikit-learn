.. include:: _contributors.rst

.. currentmodule:: sklearn

.. _changes_0_13_1:

Version 0.13.1
==============

**February 23, 2013**

The 0.13.1 release only fixes some bugs and does not add any new functionality.

Changelog
---------

- Fixed a testing error caused by the function :func:`cross_validation.train_test_split` being
  interpreted as a test by `Yaroslav Halchenko`_.

- Fixed a bug in the reassignment of small clusters in the :class:`cluster.MiniBatchKMeans`
  by `Gael Varoquaux`_.

- Fixed default value of ``gamma`` in :class:`decomposition.KernelPCA` by `Lars Buitinck`_.

- Updated joblib to ``0.7.0d`` by `Gael Varoquaux`_.

- Fixed scaling of the deviance in :class:`ensemble.GradientBoostingClassifier` by `Peter Prettenhofer`_.

- Better tie-breaking in :class:`multiclass.OneVsOneClassifier` by `Andreas Müller`_.

- Other small improvements to tests and documentation.

People
------
List of contributors for release 0.13.1 by number of commits.
 * 16  `Lars Buitinck`_
 * 12  `Andreas Müller`_
 *  8  `Gael Varoquaux`_
 *  5  Robert Marchman
 *  3  `Peter Prettenhofer`_
 *  2  Hrishikesh Huilgolkar
 *  1  Bastiaan van den Berg
 *  1  Diego Molla
 *  1  `Gilles Louppe`_
 *  1  `Mathieu Blondel`_
 *  1  `Nelle Varoquaux`_
 *  1  Rafael Cunha de Almeida
 *  1  Rolando Espinoza La fuente
 *  1  `Vlad Niculae`_
 *  1  `Yaroslav Halchenko`_


.. _changes_0_13:

Version 0.13
============

**January 21, 2013**

New Estimator Classes
---------------------

- :class:`dummy.DummyClassifier` and :class:`dummy.DummyRegressor`, two
  data-independent predictors by `Mathieu Blondel`_. Useful to sanity-check
  your estimators. See :ref:`dummy_estimators` in the user guide.
  Multioutput support added by `Arnaud Joly`_.

- :class:`decomposition.FactorAnalysis`, a transformer implementing the
  classical factor analysis, by `Christian Osendorfer`_ and `Alexandre
  Gramfort`_. See :ref:`FA` in the user guide.

- :class:`feature_extraction.FeatureHasher`, a transformer implementing the
  "hashing trick" for fast, low-memory feature extraction from string fields
  by `Lars Buitinck`_ and :class:`feature_extraction.text.HashingVectorizer`
  for text documents by `Olivier Grisel`_  See :ref:`feature_hashing` and
  :ref:`hashing_vectorizer` for the documentation and sample usage.

- :class:`pipeline.FeatureUnion`, a transformer that concatenates
  results of several other transformers by `Andreas Müller`_. See
  :ref:`feature_union` in the user guide.

- :class:`random_projection.GaussianRandomProjection`,
  :class:`random_projection.SparseRandomProjection` and the function
  :func:`random_projection.johnson_lindenstrauss_min_dim`. The first two are
  transformers implementing Gaussian and sparse random projection matrix
  by `Olivier Grisel`_ and `Arnaud Joly`_.
  See :ref:`random_projection` in the user guide.

- :class:`kernel_approximation.Nystroem`, a transformer for approximating
  arbitrary kernels by `Andreas Müller`_. See
  :ref:`nystroem_kernel_approx` in the user guide.

- :class:`preprocessing.OneHotEncoder`, a transformer that computes binary
  encodings of categorical features by `Andreas Müller`_. See
  :ref:`preprocessing_categorical_features` in the user guide.

- :class:`linear_model.PassiveAggressiveClassifier` and
  :class:`linear_model.PassiveAggressiveRegressor`, predictors implementing
  an efficient stochastic optimization for linear models by `Rob Zinkov`_ and
  `Mathieu Blondel`_. See :ref:`passive_aggressive` in the user
  guide.

- :class:`ensemble.RandomTreesEmbedding`, a transformer for creating high-dimensional
  sparse representations using ensembles of totally random trees by  `Andreas Müller`_.
  See :ref:`random_trees_embedding` in the user guide.

- :class:`manifold.SpectralEmbedding` and function
  :func:`manifold.spectral_embedding`, implementing the "laplacian
  eigenmaps" transformation for non-linear dimensionality reduction by Wei
  Li. See :ref:`spectral_embedding` in the user guide.

- :class:`isotonic.IsotonicRegression` by `Fabian Pedregosa`_, `Alexandre Gramfort`_
  and `Nelle Varoquaux`_,


Changelog
---------

- :func:`metrics.zero_one_loss` (formerly ``metrics.zero_one``) now has
  option for normalized output that reports the fraction of
  misclassifications, rather than the raw number of misclassifications. By
  Kyle Beauchamp.

- :class:`tree.DecisionTreeClassifier` and all derived ensemble models now
  support sample weighting, by `Noel Dawe`_  and `Gilles Louppe`_.

- Speedup improvement when using bootstrap samples in forests of randomized
  trees, by `Peter Prettenhofer`_  and `Gilles Louppe`_.

- Partial dependence plots for :ref:`gradient_boosting` in
  :func:`ensemble.partial_dependence.partial_dependence` by `Peter
  Prettenhofer`_. See :ref:`sphx_glr_auto_examples_ensemble_plot_partial_dependence.py` for an
  example.

- The table of contents on the website has now been made expandable by
  `Jaques Grobler`_.

- :class:`feature_selection.SelectPercentile` now breaks ties
  deterministically instead of returning all equally ranked features.

- :class:`feature_selection.SelectKBest` and
  :class:`feature_selection.SelectPercentile` are more numerically stable
  since they use scores, rather than p-values, to rank results. This means
  that they might sometimes select different features than they did
  previously.

- Ridge regression and ridge classification fitting with ``sparse_cg`` solver
  no longer has quadratic memory complexity, by `Lars Buitinck`_ and
  `Fabian Pedregosa`_.

- Ridge regression and ridge classification now support a new fast solver
  called ``lsqr``, by `Mathieu Blondel`_.

- Speed up of :func:`metrics.precision_recall_curve` by Conrad Lee.

- Added support for reading/writing svmlight files with pairwise
  preference attribute (qid in svmlight file format) in
  :func:`datasets.dump_svmlight_file` and
  :func:`datasets.load_svmlight_file` by `Fabian Pedregosa`_.

- Faster and more robust :func:`metrics.confusion_matrix` and
  :ref:`clustering_evaluation` by Wei Li.

- :func:`cross_validation.cross_val_score` now works with precomputed kernels
  and affinity matrices, by `Andreas Müller`_.

- LARS algorithm made more numerically stable with heuristics to drop
  regressors too correlated as well as to stop the path when
  numerical noise becomes predominant, by `Gael Varoquaux`_.

- Faster implementation of :func:`metrics.precision_recall_curve` by
  Conrad Lee.

- New kernel :class:`metrics.chi2_kernel` by `Andreas Müller`_, often used
  in computer vision applications.

- Fix of longstanding bug in :class:`naive_bayes.BernoulliNB` fixed by
  Shaun Jackman.

- Implemented ``predict_proba`` in :class:`multiclass.OneVsRestClassifier`,
  by Andrew Winterman.

- Improve consistency in gradient boosting: estimators
  :class:`ensemble.GradientBoostingRegressor` and
  :class:`ensemble.GradientBoostingClassifier` use the estimator
  :class:`tree.DecisionTreeRegressor` instead of the
  :class:`tree._tree.Tree` data structure by `Arnaud Joly`_.

- Fixed a floating point exception in the :ref:`decision trees <tree>`
  module, by Seberg.

- Fix :func:`metrics.roc_curve` fails when y_true has only one class
  by Wei Li.

- Add the :func:`metrics.mean_absolute_error` function which computes the
  mean absolute error. The :func:`metrics.mean_squared_error`,
  :func:`metrics.mean_absolute_error` and
  :func:`metrics.r2_score` metrics support multioutput by `Arnaud Joly`_.

- Fixed ``class_weight`` support in :class:`svm.LinearSVC` and
  :class:`linear_model.LogisticRegression` by `Andreas Müller`_. The meaning
  of ``class_weight`` was reversed as erroneously higher weight meant less
  positives of a given class in earlier releases.

- Improve narrative documentation and consistency in
  :mod:`sklearn.metrics` for regression and classification metrics
  by `Arnaud Joly`_.

- Fixed a bug in :class:`sklearn.svm.SVC` when using csr-matrices with
  unsorted indices by Xinfan Meng and `Andreas Müller`_.

- :class:`MiniBatchKMeans`: Add random reassignment of cluster centers
  with little observations attached to them, by `Gael Varoquaux`_.


API changes summary
-------------------
- Renamed all occurrences of ``n_atoms`` to ``n_components`` for consistency.
  This applies to :class:`decomposition.DictionaryLearning`,
  :class:`decomposition.MiniBatchDictionaryLearning`,
  :func:`decomposition.dict_learning`, :func:`decomposition.dict_learning_online`.

- Renamed all occurrences of ``max_iters`` to ``max_iter`` for consistency.
  This applies to :class:`semi_supervised.LabelPropagation` and
  :class:`semi_supervised.label_propagation.LabelSpreading`.

- Renamed all occurrences of ``learn_rate`` to ``learning_rate`` for
  consistency in :class:`ensemble.BaseGradientBoosting` and
  :class:`ensemble.GradientBoostingRegressor`.

- The module ``sklearn.linear_model.sparse`` is gone. Sparse matrix support
  was already integrated into the "regular" linear models.

- :func:`sklearn.metrics.mean_square_error`, which incorrectly returned the
  accumulated error, was removed. Use ``mean_squared_error`` instead.

- Passing ``class_weight`` parameters to ``fit`` methods is no longer
  supported. Pass them to estimator constructors instead.

- GMMs no longer have ``decode`` and ``rvs`` methods. Use the ``score``,
  ``predict`` or ``sample`` methods instead.

- The ``solver`` fit option in Ridge regression and classification is now
  deprecated and will be removed in v0.14. Use the constructor option
  instead.

- :class:`feature_extraction.text.DictVectorizer` now returns sparse
  matrices in the CSR format, instead of COO.

- Renamed ``k`` in :class:`cross_validation.KFold` and
  :class:`cross_validation.StratifiedKFold` to ``n_folds``, renamed
  ``n_bootstraps`` to ``n_iter`` in ``cross_validation.Bootstrap``.

- Renamed all occurrences of ``n_iterations`` to ``n_iter`` for consistency.
  This applies to :class:`cross_validation.ShuffleSplit`,
  :class:`cross_validation.StratifiedShuffleSplit`,
  :func:`utils.randomized_range_finder` and :func:`utils.randomized_svd`.

- Replaced ``rho`` in :class:`linear_model.ElasticNet` and
  :class:`linear_model.SGDClassifier` by ``l1_ratio``. The ``rho`` parameter
  had different meanings; ``l1_ratio`` was introduced to avoid confusion.
  It has the same meaning as previously ``rho`` in
  :class:`linear_model.ElasticNet` and ``(1-rho)`` in
  :class:`linear_model.SGDClassifier`.

- :class:`linear_model.LassoLars` and :class:`linear_model.Lars` now
  store a list of paths in the case of multiple targets, rather than
  an array of paths.

- The attribute ``gmm`` of :class:`hmm.GMMHMM` was renamed to ``gmm_``
  to adhere more strictly with the API.

- :func:`cluster.spectral_embedding` was moved to
  :func:`manifold.spectral_embedding`.

- Renamed ``eig_tol`` in :func:`manifold.spectral_embedding`,
  :class:`cluster.SpectralClustering` to ``eigen_tol``, renamed ``mode``
  to ``eigen_solver``.

- Renamed ``mode`` in :func:`manifold.spectral_embedding` and
  :class:`cluster.SpectralClustering` to ``eigen_solver``.

- ``classes_`` and ``n_classes_`` attributes of
  :class:`tree.DecisionTreeClassifier` and all derived ensemble models are
  now flat in case of single output problems and nested in case of
  multi-output problems.

- The ``estimators_`` attribute of
  :class:`ensemble.gradient_boosting.GradientBoostingRegressor` and
  :class:`ensemble.gradient_boosting.GradientBoostingClassifier` is now an
  array of :class:'tree.DecisionTreeRegressor'.

- Renamed ``chunk_size`` to ``batch_size`` in
  :class:`decomposition.MiniBatchDictionaryLearning` and
  :class:`decomposition.MiniBatchSparsePCA` for consistency.

- :class:`svm.SVC` and :class:`svm.NuSVC` now provide a ``classes_``
  attribute and support arbitrary dtypes for labels ``y``.
  Also, the dtype returned by ``predict`` now reflects the dtype of
  ``y`` during ``fit`` (used to be ``np.float``).

- Changed default test_size in :func:`cross_validation.train_test_split`
  to None, added possibility to infer ``test_size`` from ``train_size`` in
  :class:`cross_validation.ShuffleSplit` and
  :class:`cross_validation.StratifiedShuffleSplit`.

- Renamed function :func:`sklearn.metrics.zero_one` to
  :func:`sklearn.metrics.zero_one_loss`. Be aware that the default behavior
  in :func:`sklearn.metrics.zero_one_loss` is different from
  :func:`sklearn.metrics.zero_one`: ``normalize=False`` is changed to
  ``normalize=True``.

- Renamed function :func:`metrics.zero_one_score` to
  :func:`metrics.accuracy_score`.

- :func:`datasets.make_circles` now has the same number of inner and outer points.

- In the Naive Bayes classifiers, the ``class_prior`` parameter was moved
  from ``fit`` to ``__init__``.

People
------
List of contributors for release 0.13 by number of commits.

 * 364  `Andreas Müller`_
 * 143  `Arnaud Joly`_
 * 137  `Peter Prettenhofer`_
 * 131  `Gael Varoquaux`_
 * 117  `Mathieu Blondel`_
 * 108  `Lars Buitinck`_
 * 106  Wei Li
 * 101  `Olivier Grisel`_
 *  65  `Vlad Niculae`_
 *  54  `Gilles Louppe`_
 *  40  `Jaques Grobler`_
 *  38  `Alexandre Gramfort`_
 *  30  `Rob Zinkov`_
 *  19  Aymeric Masurelle
 *  18  Andrew Winterman
 *  17  `Fabian Pedregosa`_
 *  17  Nelle Varoquaux
 *  16  `Christian Osendorfer`_
 *  14  `Daniel Nouri`_
 *  13  :user:`Virgile Fritsch <VirgileFritsch>`
 *  13  syhw
 *  12  `Satrajit Ghosh`_
 *  10  Corey Lynch
 *  10  Kyle Beauchamp
 *   9  Brian Cheung
 *   9  Immanuel Bayer
 *   9  mr.Shu
 *   8  Conrad Lee
 *   8  `James Bergstra`_
 *   7  Tadej Janež
 *   6  Brian Cajes
 *   6  `Jake Vanderplas`_
 *   6  Michael
 *   6  Noel Dawe
 *   6  Tiago Nunes
 *   6  cow
 *   5  Anze
 *   5  Shiqiao Du
 *   4  Christian Jauvin
 *   4  Jacques Kvam
 *   4  Richard T. Guy
 *   4  `Robert Layton`_
 *   3  Alexandre Abraham
 *   3  Doug Coleman
 *   3  Scott Dickerson
 *   2  ApproximateIdentity
 *   2  John Benediktsson
 *   2  Mark Veronda
 *   2  Matti Lyra
 *   2  Mikhail Korobov
 *   2  Xinfan Meng
 *   1  Alejandro Weinstein
 *   1  `Alexandre Passos`_
 *   1  Christoph Deil
 *   1  Eugene Nizhibitsky
 *   1  Kenneth C. Arnold
 *   1  Luis Pedro Coelho
 *   1  Miroslav Batchkarov
 *   1  Pavel
 *   1  Sebastian Berg
 *   1  Shaun Jackman
 *   1  Subhodeep Moitra
 *   1  bob
 *   1  dengemann
 *   1  emanuele
 *   1  x006

