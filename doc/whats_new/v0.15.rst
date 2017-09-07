.. include:: _contributors.rst

.. currentmodule:: sklearn

.. _changes_0_15_2:

Version 0.15.2
==============

**September 4, 2014**

Bug fixes
---------

- Fixed handling of the ``p`` parameter of the Minkowski distance that was
  previously ignored in nearest neighbors models. By :user:`Nikolay
  Mayorov <nmayorov>`.

- Fixed duplicated alphas in :class:`linear_model.LassoLars` with early
  stopping on 32 bit Python. By `Olivier Grisel`_ and `Fabian Pedregosa`_.

- Fixed the build under Windows when scikit-learn is built with MSVC while
  NumPy is built with MinGW. By `Olivier Grisel`_ and :user:`Federico
  Vaggi <FedericoV>`.

- Fixed an array index overflow bug in the coordinate descent solver. By
  `Gael Varoquaux`_.

- Better handling of numpy 1.9 deprecation warnings. By `Gael Varoquaux`_.

- Removed unnecessary data copy in :class:`cluster.KMeans`.
  By `Gael Varoquaux`_.

- Explicitly close open files to avoid ``ResourceWarnings`` under Python 3.
  By Calvin Giles.

- The ``transform`` of :class:`discriminant_analysis.LinearDiscriminantAnalysis`
  now projects the input on the most discriminant directions. By Martin Billinger.

- Fixed potential overflow in ``_tree.safe_realloc`` by `Lars Buitinck`_.

- Performance optimization in :class:`isotonic.IsotonicRegression`.
  By Robert Bradshaw.

- ``nose`` is non-longer a runtime dependency to import ``sklearn``, only for
  running the tests. By `Joel Nothman`_.

- Many documentation and website fixes by `Joel Nothman`_, `Lars Buitinck`_
  :user:`Matt Pico <MattpSoftware>`, and others.

.. _changes_0_15_1:

Version 0.15.1
==============

**August 1, 2014**

Bug fixes
---------

- Made :func:`cross_validation.cross_val_score` use
  :class:`cross_validation.KFold` instead of
  :class:`cross_validation.StratifiedKFold` on multi-output classification
  problems. By :user:`Nikolay Mayorov <nmayorov>`.

- Support unseen labels :class:`preprocessing.LabelBinarizer` to restore
  the default behavior of 0.14.1 for backward compatibility. By
  :user:`Hamzeh Alsalhi <hamsal>`.

- Fixed the :class:`cluster.KMeans` stopping criterion that prevented early
  convergence detection. By Edward Raff and `Gael Varoquaux`_.

- Fixed the behavior of :class:`multiclass.OneVsOneClassifier`.
  in case of ties at the per-class vote level by computing the correct
  per-class sum of prediction scores. By `Andreas Müller`_.

- Made :func:`cross_validation.cross_val_score` and
  :class:`grid_search.GridSearchCV` accept Python lists as input data.
  This is especially useful for cross-validation and model selection of
  text processing pipelines. By `Andreas Müller`_.

- Fixed data input checks of most estimators to accept input data that
  implements the NumPy ``__array__`` protocol. This is the case for
  for ``pandas.Series`` and ``pandas.DataFrame`` in recent versions of
  pandas. By `Gael Varoquaux`_.

- Fixed a regression for :class:`linear_model.SGDClassifier` with
  ``class_weight="auto"`` on data with non-contiguous labels. By
  `Olivier Grisel`_.


.. _changes_0_15:

Version 0.15
============

**July 15, 2014**

Highlights
-----------

- Many speed and memory improvements all across the code

- Huge speed and memory improvements to random forests (and extra
  trees) that also benefit better from parallel computing.

- Incremental fit to :class:`BernoulliRBM <neural_network.BernoulliRBM>`

- Added :class:`cluster.AgglomerativeClustering` for hierarchical
  agglomerative clustering with average linkage, complete linkage and
  ward strategies.

- Added :class:`linear_model.RANSACRegressor` for robust regression
  models.

- Added dimensionality reduction with :class:`manifold.TSNE` which can be
  used to visualize high-dimensional data.


Changelog
---------

New features
............

- Added :class:`ensemble.BaggingClassifier` and
  :class:`ensemble.BaggingRegressor` meta-estimators for ensembling
  any kind of base estimator. See the :ref:`Bagging <bagging>` section of
  the user guide for details and examples. By `Gilles Louppe`_.

- New unsupervised feature selection algorithm
  :class:`feature_selection.VarianceThreshold`, by `Lars Buitinck`_.

- Added :class:`linear_model.RANSACRegressor` meta-estimator for the robust
  fitting of regression models. By :user:`Johannes Schönberger <ahojnnes>`.

- Added :class:`cluster.AgglomerativeClustering` for hierarchical
  agglomerative clustering with average linkage, complete linkage and
  ward strategies, by  `Nelle Varoquaux`_ and `Gael Varoquaux`_.

- Shorthand constructors :func:`pipeline.make_pipeline` and
  :func:`pipeline.make_union` were added by `Lars Buitinck`_.

- Shuffle option for :class:`cross_validation.StratifiedKFold`.
  By :user:`Jeffrey Blackburne <jblackburne>`.

- Incremental learning (``partial_fit``) for Gaussian Naive Bayes by
  Imran Haque.

- Added ``partial_fit`` to :class:`BernoulliRBM
  <neural_network.BernoulliRBM>`
  By :user:`Danny Sullivan <dsullivan7>`.

- Added :func:`learning_curve <learning_curve.learning_curve>` utility to
  chart performance with respect to training size. See
  :ref:`sphx_glr_auto_examples_model_selection_plot_learning_curve.py`. By Alexander Fabisch.

- Add positive option in :class:`LassoCV <linear_model.LassoCV>` and
  :class:`ElasticNetCV <linear_model.ElasticNetCV>`.
  By Brian Wignall and `Alexandre Gramfort`_.

- Added :class:`linear_model.MultiTaskElasticNetCV` and
  :class:`linear_model.MultiTaskLassoCV`. By `Manoj Kumar`_.

- Added :class:`manifold.TSNE`. By Alexander Fabisch.

Enhancements
............

- Add sparse input support to :class:`ensemble.AdaBoostClassifier` and
  :class:`ensemble.AdaBoostRegressor` meta-estimators.
  By :user:`Hamzeh Alsalhi <hamsal>`.

- Memory improvements of decision trees, by `Arnaud Joly`_.

- Decision trees can now be built in best-first manner by using ``max_leaf_nodes``
  as the stopping criteria. Refactored the tree code to use either a
  stack or a priority queue for tree building.
  By `Peter Prettenhofer`_ and `Gilles Louppe`_.

- Decision trees can now be fitted on fortran- and c-style arrays, and
  non-continuous arrays without the need to make a copy.
  If the input array has a different dtype than ``np.float32``, a fortran-
  style copy will be made since fortran-style memory layout has speed
  advantages. By `Peter Prettenhofer`_ and `Gilles Louppe`_.

- Speed improvement of regression trees by optimizing the
  the computation of the mean square error criterion. This lead
  to speed improvement of the tree, forest and gradient boosting tree
  modules. By `Arnaud Joly`_

- The ``img_to_graph`` and ``grid_tograph`` functions in
  :mod:`sklearn.feature_extraction.image` now return ``np.ndarray``
  instead of ``np.matrix`` when ``return_as=np.ndarray``.  See the
  Notes section for more information on compatibility.

- Changed the internal storage of decision trees to use a struct array.
  This fixed some small bugs, while improving code and providing a small
  speed gain. By `Joel Nothman`_.

- Reduce memory usage and overhead when fitting and predicting with forests
  of randomized trees in parallel with ``n_jobs != 1`` by leveraging new
  threading backend of joblib 0.8 and releasing the GIL in the tree fitting
  Cython code.  By `Olivier Grisel`_ and `Gilles Louppe`_.

- Speed improvement of the :mod:`sklearn.ensemble.gradient_boosting` module.
  By `Gilles Louppe`_ and `Peter Prettenhofer`_.

- Various enhancements to the  :mod:`sklearn.ensemble.gradient_boosting`
  module: a ``warm_start`` argument to fit additional trees,
  a ``max_leaf_nodes`` argument to fit GBM style trees,
  a ``monitor`` fit argument to inspect the estimator during training, and
  refactoring of the verbose code. By `Peter Prettenhofer`_.

- Faster :class:`sklearn.ensemble.ExtraTrees` by caching feature values.
  By `Arnaud Joly`_.

- Faster depth-based tree building algorithm such as decision tree,
  random forest, extra trees or gradient tree boosting (with depth based
  growing strategy) by avoiding trying to split on found constant features
  in the sample subset. By `Arnaud Joly`_.

- Add ``min_weight_fraction_leaf`` pre-pruning parameter to tree-based
  methods: the minimum weighted fraction of the input samples required to be
  at a leaf node. By `Noel Dawe`_.

- Added :func:`metrics.pairwise_distances_argmin_min`, by Philippe Gervais.

- Added predict method to :class:`cluster.AffinityPropagation` and
  :class:`cluster.MeanShift`, by `Mathieu Blondel`_.

- Vector and matrix multiplications have been optimised throughout the
  library by `Denis Engemann`_, and `Alexandre Gramfort`_.
  In particular, they should take less memory with older NumPy versions
  (prior to 1.7.2).

- Precision-recall and ROC examples now use train_test_split, and have more
  explanation of why these metrics are useful. By `Kyle Kastner`_

- The training algorithm for :class:`decomposition.NMF` is faster for
  sparse matrices and has much lower memory complexity, meaning it will
  scale up gracefully to large datasets. By `Lars Buitinck`_.

- Added svd_method option with default value to "randomized" to
  :class:`decomposition.FactorAnalysis` to save memory and
  significantly speedup computation by `Denis Engemann`_, and
  `Alexandre Gramfort`_.

- Changed :class:`cross_validation.StratifiedKFold` to try and
  preserve as much of the original ordering of samples as possible so as
  not to hide overfitting on datasets with a non-negligible level of
  samples dependency.
  By `Daniel Nouri`_ and `Olivier Grisel`_.

- Add multi-output support to :class:`gaussian_process.GaussianProcess`
  by John Novak.

- Support for precomputed distance matrices in nearest neighbor estimators
  by `Robert Layton`_ and `Joel Nothman`_.

- Norm computations optimized for NumPy 1.6 and later versions by
  `Lars Buitinck`_. In particular, the k-means algorithm no longer
  needs a temporary data structure the size of its input.

- :class:`dummy.DummyClassifier` can now be used to predict a constant
  output value. By `Manoj Kumar`_.

- :class:`dummy.DummyRegressor` has now a strategy parameter which allows
  to predict the mean, the median of the training set or a constant
  output value. By :user:`Maheshakya Wijewardena <maheshakya>`.

- Multi-label classification output in multilabel indicator format
  is now supported by :func:`metrics.roc_auc_score` and
  :func:`metrics.average_precision_score` by `Arnaud Joly`_.

- Significant performance improvements (more than 100x speedup for
  large problems) in :class:`isotonic.IsotonicRegression` by
  `Andrew Tulloch`_.

- Speed and memory usage improvements to the SGD algorithm for linear
  models: it now uses threads, not separate processes, when ``n_jobs>1``.
  By `Lars Buitinck`_.

- Grid search and cross validation allow NaNs in the input arrays so that
  preprocessors such as :class:`preprocessing.Imputer
  <preprocessing.Imputer>` can be trained within the cross validation loop,
  avoiding potentially skewed results.

- Ridge regression can now deal with sample weights in feature space
  (only sample space until then). By :user:`Michael Eickenberg <eickenberg>`.
  Both solutions are provided by the Cholesky solver.

- Several classification and regression metrics now support weighted
  samples with the new ``sample_weight`` argument:
  :func:`metrics.accuracy_score`,
  :func:`metrics.zero_one_loss`,
  :func:`metrics.precision_score`,
  :func:`metrics.average_precision_score`,
  :func:`metrics.f1_score`,
  :func:`metrics.fbeta_score`,
  :func:`metrics.recall_score`,
  :func:`metrics.roc_auc_score`,
  :func:`metrics.explained_variance_score`,
  :func:`metrics.mean_squared_error`,
  :func:`metrics.mean_absolute_error`,
  :func:`metrics.r2_score`.
  By `Noel Dawe`_.

- Speed up of the sample generator
  :func:`datasets.make_multilabel_classification`. By `Joel Nothman`_.

Documentation improvements
...........................

- The :ref:`Working With Text Data <text_data_tutorial>` tutorial
  has now been worked in to the main documentation's tutorial section.
  Includes exercises and skeletons for tutorial presentation.
  Original tutorial created by several authors including
  `Olivier Grisel`_, Lars Buitinck and many others.
  Tutorial integration into the scikit-learn documentation
  by `Jaques Grobler`_

- Added :ref:`Computational Performance <computational_performance>`
  documentation. Discussion and examples of prediction latency / throughput
  and different factors that have influence over speed. Additional tips for
  building faster models and choosing a relevant compromise between speed
  and predictive power.
  By :user:`Eustache Diemert <oddskool>`.

Bug fixes
.........

- Fixed bug in :class:`decomposition.MiniBatchDictionaryLearning` :
  ``partial_fit`` was not working properly.

- Fixed bug in :class:`linear_model.stochastic_gradient` :
  ``l1_ratio`` was used as ``(1.0 - l1_ratio)`` .

- Fixed bug in :class:`multiclass.OneVsOneClassifier` with string
  labels

- Fixed a bug in :class:`LassoCV <linear_model.LassoCV>` and
  :class:`ElasticNetCV <linear_model.ElasticNetCV>`: they would not
  pre-compute the Gram matrix with ``precompute=True`` or
  ``precompute="auto"`` and ``n_samples > n_features``. By `Manoj Kumar`_.

- Fixed incorrect estimation of the degrees of freedom in
  :func:`feature_selection.f_regression` when variates are not centered.
  By :user:`Virgile Fritsch <VirgileFritsch>`.

- Fixed a race condition in parallel processing with
  ``pre_dispatch != "all"`` (for instance, in ``cross_val_score``).
  By `Olivier Grisel`_.

- Raise error in :class:`cluster.FeatureAgglomeration` and
  :class:`cluster.WardAgglomeration` when no samples are given,
  rather than returning meaningless clustering.

- Fixed bug in :class:`gradient_boosting.GradientBoostingRegressor` with
  ``loss='huber'``: ``gamma`` might have not been initialized.

- Fixed feature importances as computed with a forest of randomized trees
  when fit with ``sample_weight != None`` and/or with ``bootstrap=True``.
  By `Gilles Louppe`_.

API changes summary
-------------------

- :mod:`sklearn.hmm` is deprecated. Its removal is planned
  for the 0.17 release.

- Use of :class:`covariance.EllipticEnvelop` has now been removed after
  deprecation.
  Please use :class:`covariance.EllipticEnvelope` instead.

- :class:`cluster.Ward` is deprecated. Use
  :class:`cluster.AgglomerativeClustering` instead.

- :class:`cluster.WardClustering` is deprecated. Use
- :class:`cluster.AgglomerativeClustering` instead.

- :class:`cross_validation.Bootstrap` is deprecated.
  :class:`cross_validation.KFold` or
  :class:`cross_validation.ShuffleSplit` are recommended instead.

- Direct support for the sequence of sequences (or list of lists) multilabel
  format is deprecated. To convert to and from the supported binary
  indicator matrix format, use
  :class:`MultiLabelBinarizer <preprocessing.MultiLabelBinarizer>`.
  By `Joel Nothman`_.

- Add score method to :class:`PCA <decomposition.PCA>` following the model of
  probabilistic PCA and deprecate
  :class:`ProbabilisticPCA <decomposition.ProbabilisticPCA>` model whose
  score implementation is not correct. The computation now also exploits the
  matrix inversion lemma for faster computation. By `Alexandre Gramfort`_.

- The score method of :class:`FactorAnalysis <decomposition.FactorAnalysis>`
  now returns the average log-likelihood of the samples. Use score_samples
  to get log-likelihood of each sample. By `Alexandre Gramfort`_.

- Generating boolean masks (the setting ``indices=False``)
  from cross-validation generators is deprecated.
  Support for masks will be removed in 0.17.
  The generators have produced arrays of indices by default since 0.10.
  By `Joel Nothman`_.

- 1-d arrays containing strings with ``dtype=object`` (as used in Pandas)
  are now considered valid classification targets. This fixes a regression
  from version 0.13 in some classifiers. By `Joel Nothman`_.

- Fix wrong ``explained_variance_ratio_`` attribute in
  :class:`RandomizedPCA <decomposition.RandomizedPCA>`.
  By `Alexandre Gramfort`_.

- Fit alphas for each ``l1_ratio`` instead of ``mean_l1_ratio`` in
  :class:`linear_model.ElasticNetCV` and :class:`linear_model.LassoCV`.
  This changes the shape of ``alphas_`` from ``(n_alphas,)`` to
  ``(n_l1_ratio, n_alphas)`` if the ``l1_ratio`` provided is a 1-D array like
  object of length greater than one.
  By `Manoj Kumar`_.

- Fix :class:`linear_model.ElasticNetCV` and :class:`linear_model.LassoCV`
  when fitting intercept and input data is sparse. The automatic grid
  of alphas was not computed correctly and the scaling with normalize
  was wrong. By `Manoj Kumar`_.

- Fix wrong maximal number of features drawn (``max_features``) at each split
  for decision trees, random forests and gradient tree boosting.
  Previously, the count for the number of drawn features started only after
  one non constant features in the split. This bug fix will affect
  computational and generalization performance of those algorithms in the
  presence of constant features. To get back previous generalization
  performance, you should modify the value of ``max_features``.
  By `Arnaud Joly`_.

- Fix wrong maximal number of features drawn (``max_features``) at each split
  for :class:`ensemble.ExtraTreesClassifier` and
  :class:`ensemble.ExtraTreesRegressor`. Previously, only non constant
  features in the split was counted as drawn. Now constant features are
  counted as drawn. Furthermore at least one feature must be non constant
  in order to make a valid split. This bug fix will affect
  computational and generalization performance of extra trees in the
  presence of constant features. To get back previous generalization
  performance, you should modify the value of ``max_features``.
  By `Arnaud Joly`_.

- Fix :func:`utils.compute_class_weight` when ``class_weight=="auto"``.
  Previously it was broken for input of non-integer ``dtype`` and the
  weighted array that was returned was wrong. By `Manoj Kumar`_.

- Fix :class:`cross_validation.Bootstrap` to return ``ValueError``
  when ``n_train + n_test > n``. By :user:`Ronald Phlypo <rphlypo>`.


People
------

List of contributors for release 0.15 by number of commits.

* 312	Olivier Grisel
* 275	Lars Buitinck
* 221	Gael Varoquaux
* 148	Arnaud Joly
* 134	Johannes Schönberger
* 119	Gilles Louppe
* 113	Joel Nothman
* 111	Alexandre Gramfort
*  95	Jaques Grobler
*  89	Denis Engemann
*  83	Peter Prettenhofer
*  83	Alexander Fabisch
*  62	Mathieu Blondel
*  60	Eustache Diemert
*  60	Nelle Varoquaux
*  49	Michael Bommarito
*  45	Manoj-Kumar-S
*  28	Kyle Kastner
*  26	Andreas Mueller
*  22	Noel Dawe
*  21	Maheshakya Wijewardena
*  21	Brooke Osborn
*  21	Hamzeh Alsalhi
*  21	Jake VanderPlas
*  21	Philippe Gervais
*  19	Bala Subrahmanyam Varanasi
*  12	Ronald Phlypo
*  10	Mikhail Korobov
*   8	Thomas Unterthiner
*   8	Jeffrey Blackburne
*   8	eltermann
*   8	bwignall
*   7	Ankit Agrawal
*   7	CJ Carey
*   6	Daniel Nouri
*   6	Chen Liu
*   6	Michael Eickenberg
*   6	ugurthemaster
*   5	Aaron Schumacher
*   5	Baptiste Lagarde
*   5	Rajat Khanduja
*   5	Robert McGibbon
*   5	Sergio Pascual
*   4	Alexis Metaireau
*   4	Ignacio Rossi
*   4	Virgile Fritsch
*   4	Sebastian Säger
*   4	Ilambharathi Kanniah
*   4	sdenton4
*   4	Robert Layton
*   4	Alyssa
*   4	Amos Waterland
*   3	Andrew Tulloch
*   3	murad
*   3	Steven Maude
*   3	Karol Pysniak
*   3	Jacques Kvam
*   3	cgohlke
*   3	cjlin
*   3	Michael Becker
*   3	hamzeh
*   3	Eric Jacobsen
*   3	john collins
*   3	kaushik94
*   3	Erwin Marsi
*   2	csytracy
*   2	LK
*   2	Vlad Niculae
*   2	Laurent Direr
*   2	Erik Shilts
*   2	Raul Garreta
*   2	Yoshiki Vázquez Baeza
*   2	Yung Siang Liau
*   2	abhishek thakur
*   2	James Yu
*   2	Rohit Sivaprasad
*   2	Roland Szabo
*   2	amormachine
*   2	Alexis Mignon
*   2	Oscar Carlsson
*   2	Nantas Nardelli
*   2	jess010
*   2	kowalski87
*   2	Andrew Clegg
*   2	Federico Vaggi
*   2	Simon Frid
*   2	Félix-Antoine Fortin
*   1	Ralf Gommers
*   1	t-aft
*   1	Ronan Amicel
*   1	Rupesh Kumar Srivastava
*   1	Ryan Wang
*   1	Samuel Charron
*   1	Samuel St-Jean
*   1	Fabian Pedregosa
*   1	Skipper Seabold
*   1	Stefan Walk
*   1	Stefan van der Walt
*   1	Stephan Hoyer
*   1	Allen Riddell
*   1	Valentin Haenel
*   1	Vijay Ramesh
*   1	Will Myers
*   1	Yaroslav Halchenko
*   1	Yoni Ben-Meshulam
*   1	Yury V. Zaytsev
*   1	adrinjalali
*   1	ai8rahim
*   1	alemagnani
*   1	alex
*   1	benjamin wilson
*   1	chalmerlowe
*   1	dzikie drożdże
*   1	jamestwebber
*   1	matrixorz
*   1	popo
*   1	samuela
*   1	François Boulogne
*   1	Alexander Measure
*   1	Ethan White
*   1	Guilherme Trein
*   1	Hendrik Heuer
*   1	IvicaJovic
*   1	Jan Hendrik Metzen
*   1	Jean Michel Rouly
*   1	Eduardo Ariño de la Rubia
*   1	Jelle Zijlstra
*   1	Eddy L O Jansson
*   1	Denis
*   1	John
*   1	John Schmidt
*   1	Jorge Cañardo Alastuey
*   1	Joseph Perla
*   1	Joshua Vredevoogd
*   1	José Ricardo
*   1	Julien Miotte
*   1	Kemal Eren
*   1	Kenta Sato
*   1	David Cournapeau
*   1	Kyle Kelley
*   1	Daniele Medri
*   1	Laurent Luce
*   1	Laurent Pierron
*   1	Luis Pedro Coelho
*   1	DanielWeitzenfeld
*   1	Craig Thompson
*   1	Chyi-Kwei Yau
*   1	Matthew Brett
*   1	Matthias Feurer
*   1	Max Linke
*   1	Chris Filo Gorgolewski
*   1	Charles Earl
*   1	Michael Hanke
*   1	Michele Orrù
*   1	Bryan Lunt
*   1	Brian Kearns
*   1	Paul Butler
*   1	Paweł Mandera
*   1	Peter
*   1	Andrew Ash
*   1	Pietro Zambelli
*   1	staubda

