.. currentmodule:: sklearn

.. _changes_0_16:

0.16
====


Changelog
---------

New features
............

   - Incremental fit for :class:`GaussianNB <naive_bayes.GaussianNB>`.

   - Add ``sample_weight`` support to :class:`dummy.DummyClassifier`. By
     `Arnaud Joly`_.

   - Add the :func:`metrics.label_ranking_average_precision_score` metrics. By
     `Arnaud Joly`_.

   - Added :class:`linear_model.LogisticRegressionCV`. By
     `Manoj Kumar`_, `Fabian Pedregosa`_, `Gael Varoquaux`_
     and `Alexandre Gramfort`_.

   - Added ``warm_start`` constructor parameter to make it possible for any
     trained forest model to grow additional trees incrementally. By
     `Laurent Direr`_.

   - Add ``sample_weight`` support to :class:`ensemble.GradientBoostingClassifier` and
     :class:`ensemble.GradientBoostingRegressor`. By
     `Peter Prettenhofer`_.

   - Added :class:`decomposition.IncrementalPCA`, an implementation of the PCA
     algorithm that supports out-of-core learning with a ``partial_fit``
     method. By `Kyle Kastner`_.


Enhancements
............

   - Add support for sample weights in scorer objects.  Metrics with sample
     weight support will automatically benefit from it. By `Noel Dawe`_ and
     `Vlad Niculae`_.

   - Added ``newton-cg`` and `lbfgs` solver support in
     :class:`linear_model.LogisticRegression`. By `Manoj Kumar`_.

   - Add ``selection="random"`` parameter to implement stochastic coordinate
     descent for :class:`linear_model.Lasso`, :class:`linear_model.ElasticNet`
     and related. By `Manoj Kumar`_.

   - Add ``sample_weight`` parameter to `metrics.jaccard_similarity_score` and
     `metrics.log_loss`. By `Jatin Shah`_.

   - Support sparse multilabel indicator representation in
     :class:`preprocessing.LabelBinarizer` and
     :class:`multiclass.OneVsRestClassifier` (by `Hamzeh Alsalhi`_ with thanks
     to `Rohit Sivaprasad`_), as well as evaluation metrics (by
     `Joel Nothman`_).

   - Add ``multi_class="multinomial"`` option in
     :class:`linear_model.LogisticRegression` to implement a Logistic
     Regression solver that minimizes the cross-entropy or multinomial loss
     instead of the default One-vs-Rest setting. By `Lars Buitinck`_ and
     `Manoj Kumar`_.

   - ``DictVectorizer`` can now perform ``fit_transform`` on an iterable in a
     single pass, when giving the option ``sort=False``. By Dan Blanchard.

   - :class:`GridSearchCV` and :class:`RandomizedSearchCV` can now be
     configured to work with estimators that may fail and raise errors on
     individual folds. This option is controlled by the `error_score`
     parameter. This does not affect errors raised on re-fit. By 
	 `Michal Romaniuk`_.


Documentation improvements
..........................


Bug fixes
.........

    - The :class:`decomposition.PCA` now undoes whitening in its
     ``inverse_transform``. Also, its ``components_`` now always have unit
     length. By Michael Eickenberg.

    - Fix incomplete download of the dataset when
      :func:`datasets.download_20newsgroups` is called. By `Manoj Kumar`_.

    - Various fixes to the Gaussian processes subpackage by Vincent Dubourg
      and Jan Hendrik Metzen.

    - Calling ``partial_fit`` with ``class_weight=='auto'`` throws an
      appropriate error message and suggests a work around.
      By `Danny Sullivan`_.

    - :class:`RBFSampler <kernel_approximation.RBFSampler>` with ``gamma=g``
      formerly approximated :func:`rbf_kernel <metrics.pairwise.rbf_kernel>`
      with ``gamma=g/2.``; the definition of ``gamma`` is now consistent,
      which may substantially change your results if you use a fixed value.
      (If you cross-validated over ``gamma``, it probably doesn't matter
      too much.) By `Dougal Sutherland`_.


API changes summary
-------------------

    - :class:`GridSearchCV <grid_search.GridSearchCV>` and
      :func:`cross_val_score <cross_validation.cross_val_score>` and other
      meta-estimators don't convert pandas DataFrames into arrays any more,
      allowing DataFrame specific operations in custom estimators.

    - :func:`multiclass.fit_ovr`, :func:`multiclass.predict_ovr`,
      :func:`predict_proba_ovr`,
      :func:`multiclass.fit_ovo`, :func:`multiclass.predict_ovo`,
      :func:`multiclass.fit_ecoc` and :func:`multiclass.predict_ecoc`
      are deprecated. Use the underlying estimators instead.

    - Nearest neighbors estimators used to take arbitrary keyword arguments
      and pass these to their distance metric. This will no longer be supported
      in scikit-learn 0.18; use the ``metric_params`` argument instead.

    - `n_jobs` parameter of the fit method shifted to the constructor of the
       LinearRegression class.

    - The ``predict_proba`` method of :class:`multiclass.OneVsRestClassifier`
      now returns two probabilities per sample in the multiclass case; this
      is consistent with other estimators and with the method's documentation,
      but previous versions accidentally returned only the positive
      probability. Fixed by Will Lamond and `Lars Buitinck`_.


.. _changes_0_15_2:

0.15.2
======

Bug fixes
---------

  - Fixed handling of the ``p`` parameter of the Minkowski distance that was
    previously ignored in nearest neighbors models. By `Nikolay Mayorov`_.

  - Fixed duplicated alphas in :class:`linear_model.LassoLars` with early
    stopping on 32 bit Python. By `Olivier Grisel`_ and `Fabian Pedregosa`_.

  - Fixed the build under Windows when scikit-learn is built with MSVC while
    NumPy is built with MinGW. By `Olivier Grisel`_ and Federico Vaggi.

  - Fixed an array index overflow bug in the coordinate descent solver. By
    `Gael Varoquaux`_.

  - Better handling of numpy 1.9 deprecation warnings. By `Gael Varoquaux`_.

  - Removed unnecessary data copy in :class:`cluster.KMeans`.
    By `Gael Varoquaux`_.

  - Explicitly close open files to avoid ``ResourceWarnings`` under Python 3.
    By Calvin Giles.

  - The ``transform`` of :class:`lda.LDA` now projects the input on the most
    discriminant directions. By Martin Billinger.

  - Fixed potential overflow in ``_tree.safe_realloc`` by `Lars Buitinck`_.

  - Performance optimization in :class:`isotonic.IsotonicRegression`.
    By Robert Bradshaw.

  - ``nose`` is non-longer a runtime dependency to import ``sklearn``, only for
    running the tests. By `Joel Nothman`_.

  - Many documentation and website fixes by `Joel Nothman`_, `Lars Buitinck`_
    and others.

.. _changes_0_15_1:

0.15.1
======

Bug fixes
---------

   - Made :func:`cross_validation.cross_val_score` use
     :class:`cross_validation.KFold` instead of
     :class:`cross_validation.StratifiedKFold` on multi-output classification
     problems. By `Nikolay Mayorov`_.

   - Support unseen labels :class:`preprocessing.LabelBinarizer` to restore
     the default behavior of 0.14.1 for backward compatibility. By
     `Hamzeh Alsalhi`_.

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

0.15
====

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
     fitting of regression models. By Johannes Schönberger.

   - Added :class:`cluster.AgglomerativeClustering` for hierarchical
     agglomerative clustering with average linkage, complete linkage and
     ward strategies, by  `Nelle Varoquaux`_ and `Gael Varoquaux`_.

   - Shorthand constructors :func:`pipeline.make_pipeline` and
     :func:`pipeline.make_union` were added by `Lars Buitinck`_.

   - Shuffle option for :class:`cross_validation.StratifiedKFold`.
     By `Jeffrey Blackburne`_.

   - Incremental learning (``partial_fit``) for Gaussian Naive Bayes by
     Imran Haque.

   - Added ``partial_fit`` to :class:`BernoulliRBM
     <neural_network.BernoulliRBM>`
     By `Danny Sullivan`_.

   - Added :func:`learning_curve <learning_curve.learning_curve>` utility to
     chart performance with respect to training size. See
     :ref:`example_model_selection_plot_learning_curve.py`. By Alexander Fabisch.

   - Add positive option in :class:`LassoCV <linear_model.LassoCV>` and
     :class:`ElasticNetCV <linear_model.ElasticNetCV>`.
     By Brian Wignall and `Alexandre Gramfort`_.

   - Added :class:`linear_model.MultiTaskElasticNetCV` and
     :class:`linear_model.MultiTaskLassoCV`. By `Manoj Kumar`_.

Enhancements
............

   - Add sparse input support to :class:`ensemble.AdaBoostClassifier` and
     :class:`ensemble.AdaBoostRegressor` meta-estimators.
     By `Hamzeh Alsalhi`_.

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

   - Norm computations optimized for NumPy 1.6 and later versions by
     `Lars Buitinck`_. In particular, the k-means algorithm no longer
     needs a temporary data structure the size of its input.

   - :class:`dummy.DummyClassifier` can now be used to predict a constant
     output value. By `Manoj Kumar`_.

   - :class:`dummy.DummyRegressor` has now a strategy parameter which allows
     to predict the mean, the median of the training set or a constant
     output value. By `Maheshakya Wijewardena`_.

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
     (only sample space until then). By `Michael Eickenberg`_.
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
     By `Eustache Diemert`_.

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
     By `Virgile Fritsch`_.

   - Fixed a race condition in parallel processing with
     ``pre_dispatch != "all"`` (for instance in ``cross_val_score``).
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
     when ``n_train + n_test > n``. By `Ronald Phlypo`_.


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
*   4	Sebastian Saeger
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


.. _changes_0_14:

0.14
=======

Changelog
---------

   - Missing values with sparse and dense matrices can be imputed with the
     transformer :class:`preprocessing.Imputer` by `Nicolas Trésegnie`_.

   - The core implementation of decisions trees has been rewritten from
     scratch, allowing for faster tree induction and lower memory
     consumption in all tree-based estimators. By `Gilles Louppe`_.

   - Added :class:`ensemble.AdaBoostClassifier` and
     :class:`ensemble.AdaBoostRegressor`, by `Noel Dawe`_  and
     `Gilles Louppe`_. See the :ref:`AdaBoost <adaboost>` section of the user
     guide for details and examples.

   - Added :class:`grid_search.RandomizedSearchCV` and
     :class:`grid_search.ParameterSampler` for randomized hyperparameter
     optimization. By `Andreas Müller`_.

   - Added :ref:`biclustering <biclustering>` algorithms
     (:class:`sklearn.cluster.bicluster.SpectralCoclustering` and
     :class:`sklearn.cluster.bicluster.SpectralBiclustering`), data
     generation methods (:func:`sklearn.datasets.make_biclusters` and
     :func:`sklearn.datasets.make_checkerboard`), and scoring metrics
     (:func:`sklearn.metrics.consensus_score`). By `Kemal Eren`_.

   - Added :ref:`Restricted Boltzmann Machines<rbm>`
     (:class:`neural_network.BernoulliRBM`). By `Yann Dauphin`_.

   - Python 3 support by `Justin Vincent`_, `Lars Buitinck`_,
     `Subhodeep Moitra`_ and `Olivier Grisel`_. All tests now pass under
     Python 3.3.

   - Ability to pass one penalty (alpha value) per target in
     :class:`linear_model.Ridge`, by @eickenberg and `Mathieu Blondel`_.

   - Fixed :mod:`sklearn.linear_model.stochastic_gradient.py` L2 regularization
     issue (minor practical significance).
     By `Norbert Crombach`_ and `Mathieu Blondel`_ .

   - Added an interactive version of `Andreas Müller`_'s
     `Machine Learning Cheat Sheet (for scikit-learn)
     <http://peekaboo-vision.blogspot.de/2013/01/machine-learning-cheat-sheet-for-scikit.html>`_
     to the documentation. See :ref:`Choosing the right estimator <ml_map>`.
     By `Jaques Grobler`_.

   - :class:`grid_search.GridSearchCV` and
     :func:`cross_validation.cross_val_score` now support the use of advanced
     scoring function such as area under the ROC curve and f-beta scores.
     See :ref:`scoring_parameter` for details. By `Andreas Müller`_
     and `Lars Buitinck`_.
     Passing a function from :mod:`sklearn.metrics` as ``score_func`` is
     deprecated.

   - Multi-label classification output is now supported by
     :func:`metrics.accuracy_score`, :func:`metrics.zero_one_loss`,
     :func:`metrics.f1_score`, :func:`metrics.fbeta_score`,
     :func:`metrics.classification_report`,
     :func:`metrics.precision_score` and :func:`metrics.recall_score`
     by `Arnaud Joly`_.

   - Two new metrics :func:`metrics.hamming_loss` and
     :func:`metrics.jaccard_similarity_score`
     are added with multi-label support by `Arnaud Joly`_.

   - Speed and memory usage improvements in
     :class:`feature_extraction.text.CountVectorizer` and
     :class:`feature_extraction.text.TfidfVectorizer`,
     by Jochen Wersdörfer and Roman Sinayev.

   - The ``min_df`` parameter in
     :class:`feature_extraction.text.CountVectorizer` and
     :class:`feature_extraction.text.TfidfVectorizer`, which used to be 2,
     has been reset to 1 to avoid unpleasant surprises (empty vocabularies)
     for novice users who try it out on tiny document collections.
     A value of at least 2 is still recommended for practical use.

   - :class:`svm.LinearSVC`, :class:`linear_model.SGDClassifier` and
     :class:`linear_model.SGDRegressor` now have a ``sparsify`` method that
     converts their ``coef_`` into a sparse matrix, meaning stored models
     trained using these estimators can be made much more compact.

   - :class:`linear_model.SGDClassifier` now produces multiclass probability
     estimates when trained under log loss or modified Huber loss.

   - Hyperlinks to documentation in example code on the website by
     `Martin Luessi`_.

   - Fixed bug in :class:`preprocessing.MinMaxScaler` causing incorrect scaling
     of the features for non-default ``feature_range`` settings. By `Andreas
     Müller`_.

   - ``max_features`` in :class:`tree.DecisionTreeClassifier`,
     :class:`tree.DecisionTreeRegressor` and all derived ensemble estimators
     now supports percentage values. By `Gilles Louppe`_.

   - Performance improvements in :class:`isotonic.IsotonicRegression` by
     `Nelle Varoquaux`_.

   - :func:`metrics.accuracy_score` has an option normalize to return
     the fraction or the number of correctly classified sample
     by `Arnaud Joly`_.

   - Added :func:`metrics.log_loss` that computes log loss, aka cross-entropy
     loss. By Jochen Wersdörfer and `Lars Buitinck`_.

   - A bug that caused :class:`ensemble.AdaBoostClassifier`'s to output
     incorrect probabilities has been fixed.

   - Feature selectors now share a mixin providing consistent ``transform``,
     ``inverse_transform`` and ``get_support`` methods. By `Joel Nothman`_.

   - A fitted :class:`grid_search.GridSearchCV` or
     :class:`grid_search.RandomizedSearchCV` can now generally be pickled.
     By `Joel Nothman`_.

   - Refactored and vectorized implementation of :func:`metrics.roc_curve`
     and :func:`metrics.precision_recall_curve`. By `Joel Nothman`_.

   - The new estimator :class:`sklearn.decomposition.TruncatedSVD`
     performs dimensionality reduction using SVD on sparse matrices,
     and can be used for latent semantic analysis (LSA).
     By `Lars Buitinck`_.

   - Added self-contained example of out-of-core learning on text data
     :ref:`example_applications_plot_out_of_core_classification.py`.
     By `Eustache Diemert`_.

   - The default number of components for
     :class:`sklearn.decomposition.RandomizedPCA` is now correctly documented
     to be ``n_features``. This was the default behavior, so programs using it
     will continue to work as they did.

   - :class:`sklearn.cluster.KMeans` now fits several orders of magnitude
     faster on sparse data (the speedup depends on the sparsity). By
     `Lars Buitinck`_.

   - Reduce memory footprint of FastICA by `Denis Engemann`_ and
     `Alexandre Gramfort`_.

   - Verbose output in :mod:`sklearn.ensemble.gradient_boosting` now uses
     a column format and prints progress in decreasing frequency.
     It also shows the remaining time. By `Peter Prettenhofer`_.

   - :mod:`sklearn.ensemble.gradient_boosting` provides out-of-bag improvement
     :attr:`~sklearn.ensemble.GradientBoostingRegressor.oob_improvement_`
     rather than the OOB score for model selection. An example that shows
     how to use OOB estimates to select the number of trees was added.
     By `Peter Prettenhofer`_.

   - Most metrics now support string labels for multiclass classification
     by `Arnaud Joly`_ and `Lars Buitinck`_.

   - New OrthogonalMatchingPursuitCV class by `Alexandre Gramfort`_
     and `Vlad Niculae`_.

   - Fixed a bug in :class:`sklearn.covariance.GraphLassoCV`: the
     'alphas' parameter now works as expected when given a list of
     values. By Philippe Gervais.

   - Fixed an important bug in :class:`sklearn.covariance.GraphLassoCV`
     that prevented all folds provided by a CV object to be used (only
     the first 3 were used). When providing a CV object, execution
     time may thus increase significantly compared to the previous
     version (bug results are correct now). By Philippe Gervais.

   - :class:`cross_validation.cross_val_score` and the :mod:`grid_search`
     module is now tested with multi-output data by `Arnaud Joly`_.

   - :func:`datasets.make_multilabel_classification` can now return
     the output in label indicator multilabel format  by `Arnaud Joly`_.

   - K-nearest neighbors, :class:`neighbors.KNeighborsRegressor`
     and :class:`neighbors.RadiusNeighborsRegressor`,
     and radius neighbors, :class:`neighbors.RadiusNeighborsRegressor` and
     :class:`neighbors.RadiusNeighborsClassifier` support multioutput data
     by `Arnaud Joly`_.

   - Random state in LibSVM-based estimators (:class:`svm.SVC`, :class:`NuSVC`,
     :class:`OneClassSVM`, :class:`svm.SVR`, :class:`svm.NuSVR`) can now be
     controlled.  This is useful to ensure consistency in the probability
     estimates for the classifiers trained with ``probability=True``. By
     `Vlad Niculae`_.

   - Out-of-core learning support for discrete naive Bayes classifiers
     :class:`sklearn.naive_bayes.MultinomialNB` and
     :class:`sklearn.naive_bayes.BernoulliNB` by adding the ``partial_fit``
     method by `Olivier Grisel`_.

   - New website design and navigation by `Gilles Louppe`_, `Nelle Varoquaux`_,
     Vincent Michel and `Andreas Müller`_.

   - Improved documentation on :ref:`multi-class, multi-label and multi-output
     classification <multiclass>` by `Yannick Schwartz`_ and `Arnaud Joly`_.

   - Better input and error handling in the :mod:`metrics` module by
     `Arnaud Joly`_ and `Joel Nothman`_.

   - Speed optimization of the :mod:`hmm` module by `Mikhail Korobov`_

   - Significant speed improvements for :class:`sklearn.cluster.DBSCAN`
     by `cleverless <https://github.com/cleverless>`_


API changes summary
-------------------

   - The :func:`auc_score` was renamed :func:`roc_auc_score`.

   - Testing scikit-learn with ``sklearn.test()`` is deprecated. Use
     ``nosetests sklearn`` from the command line.

   - Feature importances in :class:`tree.DecisionTreeClassifier`,
     :class:`tree.DecisionTreeRegressor` and all derived ensemble estimators
     are now computed on the fly when accessing  the ``feature_importances_``
     attribute. Setting ``compute_importances=True`` is no longer required.
     By `Gilles Louppe`_.

   - :class:`linear_model.lasso_path` and
     :class:`linear_model.enet_path` can return its results in the same
     format as that of :class:`linear_model.lars_path`. This is done by
     setting the ``return_models`` parameter to ``False``. By
     `Jaques Grobler`_ and `Alexandre Gramfort`_

   - :class:`grid_search.IterGrid` was renamed to
     :class:`grid_search.ParameterGrid`.

   - Fixed bug in :class:`KFold` causing imperfect class balance in some
     cases. By `Alexandre Gramfort`_ and Tadej Janež.

   - :class:`sklearn.neighbors.BallTree` has been refactored, and a
     :class:`sklearn.neighbors.KDTree` has been
     added which shares the same interface.  The Ball Tree now works with
     a wide variety of distance metrics.  Both classes have many new
     methods, including single-tree and dual-tree queries, breadth-first
     and depth-first searching, and more advanced queries such as
     kernel density estimation and 2-point correlation functions.
     By `Jake Vanderplas`_

   - Support for scipy.spatial.cKDTree within neighbors queries has been
     removed, and the functionality replaced with the new :class:`KDTree`
     class.

   - :class:`sklearn.neighbors.KernelDensity` has been added, which performs
     efficient kernel density estimation with a variety of kernels.

   - :class:`sklearn.decomposition.KernelPCA` now always returns output with
     ``n_components`` components, unless the new parameter ``remove_zero_eig``
     is set to ``True``. This new behavior is consistent with the way
     kernel PCA was always documented; previously, the removal of components
     with zero eigenvalues was tacitly performed on all data.

   - ``gcv_mode="auto"`` no longer tries to perform SVD on a densified
     sparse matrix in :class:`sklearn.linear_model.RidgeCV`.

   - Sparse matrix support in :class:`sklearn.decomposition.RandomizedPCA`
     is now deprecated in favor of the new ``TruncatedSVD``.

   - :class:`cross_validation.KFold` and
     :class:`cross_validation.StratifiedKFold` now enforce `n_folds >= 2`
     otherwise a ``ValueError`` is raised. By `Olivier Grisel`_.

   - :func:`datasets.load_files`'s ``charset`` and ``charset_errors``
     parameters were renamed ``encoding`` and ``decode_errors``.

   - Attribute ``oob_score_`` in :class:`sklearn.ensemble.GradientBoostingRegressor`
     and :class:`sklearn.ensemble.GradientBoostingClassifier`
     is deprecated and has been replaced by ``oob_improvement_`` .

   - Attributes in OrthogonalMatchingPursuit have been deprecated
     (copy_X, Gram, ...) and precompute_gram renamed precompute
     for consistency. See #2224.

   - :class:`sklearn.preprocessing.StandardScaler` now converts integer input
     to float, and raises a warning. Previously it rounded for dense integer
     input.

   - :class:`sklearn.multiclass.OneVsRestClassifier` now has a
     ``decision_function`` method. This will return the distance of each
     sample from the decision boundary for each class, as long as the
     underlying estimators implement the ``decision_function`` method.
     By `Kyle Kastner`_.

   - Better input validation, warning on unexpected shapes for y.

People
------
List of contributors for release 0.14 by number of commits.

 * 277  Gilles Louppe
 * 245  Lars Buitinck
 * 187  Andreas Mueller
 * 124  Arnaud Joly
 * 112  Jaques Grobler
 * 109  Gael Varoquaux
 * 107  Olivier Grisel
 * 102  Noel Dawe
 *  99  Kemal Eren
 *  79  Joel Nothman
 *  75  Jake VanderPlas
 *  73  Nelle Varoquaux
 *  71  Vlad Niculae
 *  65  Peter Prettenhofer
 *  64  Alexandre Gramfort
 *  54  Mathieu Blondel
 *  38  Nicolas Trésegnie
 *  35  eustache
 *  27  Denis Engemann
 *  25  Yann N. Dauphin
 *  19  Justin Vincent
 *  17  Robert Layton
 *  15  Doug Coleman
 *  14  Michael Eickenberg
 *  13  Robert Marchman
 *  11  Fabian Pedregosa
 *  11  Philippe Gervais
 *  10  Jim Holmström
 *  10  Tadej Janež
 *  10  syhw
 *   9  Mikhail Korobov
 *   9  Steven De Gryze
 *   8  sergeyf
 *   7  Ben Root
 *   7  Hrishikesh Huilgolkar
 *   6  Kyle Kastner
 *   6  Martin Luessi
 *   6  Rob Speer
 *   5  Federico Vaggi
 *   5  Raul Garreta
 *   5  Rob Zinkov
 *   4  Ken Geis
 *   3  A. Flaxman
 *   3  Denton Cockburn
 *   3  Dougal Sutherland
 *   3  Ian Ozsvald
 *   3  Johannes Schönberger
 *   3  Robert McGibbon
 *   3  Roman Sinayev
 *   3  Szabo Roland
 *   2  Diego Molla
 *   2  Imran Haque
 *   2  Jochen Wersdörfer
 *   2  Sergey Karayev
 *   2  Yannick Schwartz
 *   2  jamestwebber
 *   1  Abhijeet Kolhe
 *   1  Alexander Fabisch
 *   1  Bastiaan van den Berg
 *   1  Benjamin Peterson
 *   1  Daniel Velkov
 *   1  Fazlul Shahriar
 *   1  Felix Brockherde
 *   1  Félix-Antoine Fortin
 *   1  Harikrishnan S
 *   1  Jack Hale
 *   1  JakeMick
 *   1  James McDermott
 *   1  John Benediktsson
 *   1  John Zwinck
 *   1  Joshua Vredevoogd
 *   1  Justin Pati
 *   1  Kevin Hughes
 *   1  Kyle Kelley
 *   1  Matthias Ekman
 *   1  Miroslav Shubernetskiy
 *   1  Naoki Orii
 *   1  Norbert Crombach
 *   1  Rafael Cunha de Almeida
 *   1  Rolando Espinoza La fuente
 *   1  Seamus Abshere
 *   1  Sergey Feldman
 *   1  Sergio Medina
 *   1  Stefano Lattarini
 *   1  Steve Koch
 *   1  Sturla Molden
 *   1  Thomas Jarosch
 *   1  Yaroslav Halchenko

.. _changes_0_13_1:

0.13.1
======

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

0.13
====

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
     Prettenhofer`_. See :ref:`example_ensemble_plot_partial_dependence.py` for an
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
 *  13  `Virgile Fritsch`_
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


.. _changes_0_12.1:

0.12.1
=======

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

 - Fix MinCovDet breaking with X.shape = (3, 1) by `Virgile Fritsch`_

 - Fix clone of SGD objects by `Peter Prettenhofer`_

 - Stabilize GMM by `Virgile Fritsch`_

People
------

 *  14  `Peter Prettenhofer`_
 *  12  `Gael Varoquaux`_
 *  10  `Andreas Müller`_
 *   5  `Lars Buitinck`_
 *   3  `Virgile Fritsch`_
 *   1  `Alexandre Gramfort`_
 *   1  `Gilles Louppe`_
 *   1  `Mathieu Blondel`_

.. _changes_0_12:

0.12
====

Changelog
---------

   - Various speed improvements of the :ref:`decision trees <tree>` module, by
     `Gilles Louppe`_.

   - :class:`ensemble.GradientBoostingRegressor` and
     :class:`ensemble.GradientBoostingClassifier` now support feature subsampling
     via the ``max_features`` argument, by `Peter Prettenhofer`_.

   - Added Huber and Quantile loss functions to
     :class:`ensemble.GradientBoostingRegressor`, by `Peter Prettenhofer`_.

   - :ref:`Decision trees <tree>` and :ref:`forests of randomized trees <forest>`
     now support multi-output classification and regression problems, by
     `Gilles Louppe`_.

   - Added :class:`preprocessing.LabelEncoder`, a simple utility class to
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

   - Added :func:`metrics.auc_score` and
     :func:`metrics.average_precision_score` convenience functions by `Andreas
     Müller`_.

   - Improved sparse matrix support in the :ref:`feature_selection`
     module by `Andreas Müller`_.

   - New word boundaries-aware character n-gram analyzer for the
     :ref:`text_feature_extraction` module by `@kernc`_.

   - Fixed bug in spectral clustering that led to single point clusters
     by `Andreas Müller`_.

   - In :class:`feature_extraction.text.CountVectorizer`, added an option to
     ignore infrequent words, ``min_df`` by  `Andreas Müller`_.

   - Add support for multiple targets in some linear models (ElasticNet, Lasso
     and OrthogonalMatchingPursuit) by `Vlad Niculae`_ and
     `Alexandre Gramfort`_.

   - Fixes in :class:`decomposition.ProbabilisticPCA` score function by Wei Li.

   - Fixed feature importance computation in
     :ref:`gradient_boosting`.

API changes summary
-------------------

   - The old ``scikits.learn`` package has disappeared; all code should import
     from ``sklearn`` instead, which was introduced in 0.9.

   - In :func:`metrics.roc_curve`, the ``thresholds`` array is now returned
     with it's order reversed, in order to keep it consistent with the order
     of the returned ``fpr`` and ``tpr``.

   - In :class:`hmm` objects, like :class:`hmm.GaussianHMM`,
     :class:`hmm.MultinomialHMM`, etc., all parameters must be passed to the
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

   - In :class:`feature_extraction.text.CountVectorizer` the parameters
     ``min_n`` and ``max_n`` were joined to the parameter ``n_gram_range`` to
     enable grid-searching both at once.

   - In :class:`feature_extraction.text.CountVectorizer`, words that appear
     only in one document are now ignored by default. To reproduce
     the previous behavior, set ``min_df=1``.

   - Fixed API inconsistency: :meth:`linear_model.SGDClassifier.predict_proba` now
     returns 2d array when fit on two classes.

   - Fixed API inconsistency: :meth:`qda.QDA.decision_function` and
     :meth:`lda.LDA.decision_function` now return 1d arrays when fit on two
     classes.

   - Grid of alphas used for fitting :class:`linear_model.LassoCV` and
     :class:`linear_model.ElasticNetCV` is now stored
     in the attribute ``alphas_`` rather than overriding the init parameter
     ``alphas``.

   - Linear models when alpha is estimated by cross-validation store
     the estimated value in the ``alpha_`` attribute rather than just
     ``alpha`` or ``best_alpha``.

   - :class:`ensemble.GradientBoostingClassifier` now supports
     :meth:`ensemble.GradientBoostingClassifier.staged_predict_proba`, and
     :meth:`ensemble.GradientBoostingClassifier.staged_predict`.

   - :class:`svm.sparse.SVC` and other sparse SVM classes are now deprecated.
     The all classes in the :ref:`svm` module now automatically select the
     sparse or dense representation base on the input.

   - All clustering algorithms now interpret the array ``X`` given to ``fit`` as
     input data, in particular :class:`cluster.SpectralClustering` and
     :class:`cluster.AffinityPropagation` which previously expected affinity matrices.

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
 *  12  `@kernc`_
 *  11  `Virgile Fritsch`_
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

0.11
====

Changelog
---------

Highlights
.............

   - Gradient boosted regression trees (:ref:`gradient_boosting`)
     for classification and regression by `Peter Prettenhofer`_
     and `Scott White`_ .

   - Simple dict-based feature loader with support for categorical variables
     (:class:`feature_extraction.DictVectorizer`) by `Lars Buitinck`_.

   - Added Matthews correlation coefficient (:func:`metrics.matthews_corrcoef`)
     and added macro and micro average options to
     :func:`metrics.precision_score`, :func:`metrics.recall_score` and
     :func:`metrics.f1_score` by `Satrajit Ghosh`_.

   - :ref:`out_of_bag` of generalization error for :ref:`ensemble`
     by `Andreas Müller`_.

   - :ref:`randomized_l1`: Randomized sparse linear models for feature
     selection, by `Alexandre Gramfort`_ and `Gael Varoquaux`_

   - :ref:`label_propagation` for semi-supervised learning, by Clay
     Woolam. **Note** the semi-supervised API is still work in progress,
     and may change.

   - Added BIC/AIC model selection to classical :ref:`gmm` and unified
     the API with the remainder of scikit-learn, by `Bertrand Thirion`_

   - Added :class:`sklearn.cross_validation.StratifiedShuffleSplit`, which is
     a :class:`sklearn.cross_validation.ShuffleSplit` with balanced splits,
     by Yannick Schwartz.

   - :class:`sklearn.neighbors.NearestCentroid` classifier added, along with a
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
     :class:`linear_model.LogisticRegression` merged by `Lars Buitinck`_.

   - Regressors can now be used as base estimator in the :ref:`multiclass`
     module by `Mathieu Blondel`_.

   - Added n_jobs option to :func:`metrics.pairwise.pairwise_distances`
     and :func:`metrics.pairwise.pairwise_kernels` for parallel computation,
     by `Mathieu Blondel`_.

   - :ref:`k_means` can now be run in parallel, using the ``n_jobs`` argument
     to either :ref:`k_means` or :class:`KMeans`, by `Robert Layton`_.

   - Improved :ref:`cross_validation` and :ref:`grid_search` documentation
     and introduced the new :func:`cross_validation.train_test_split`
     helper function by `Olivier Grisel`_

   - :class:`svm.SVC` members ``coef_`` and ``intercept_`` changed sign for
     consistency with ``decision_function``; for ``kernel==linear``,
     ``coef_`` was fixed in the the one-vs-one case, by `Andreas Müller`_.

   - Performance improvements to efficient leave-one-out cross-validated
     Ridge regression, esp. for the ``n_samples > n_features`` case, in
     :class:`linear_model.RidgeCV`, by Reuben Fletcher-Costin.

   - Refactoring and simplification of the :ref:`text_feature_extraction`
     API and fixed a bug that caused possible negative IDF,
     by `Olivier Grisel`_.

   - Beam pruning option in :class:`_BaseHMM` module has been removed since it
     is difficult to Cythonize. If you are interested in contributing a Cython
     version, you can use the python version in the git history as a reference.

   - Classes in :ref:`neighbors` now support arbitrary Minkowski metric for
     nearest neighbors searches. The metric can be specified by argument ``p``.

API changes summary
-------------------

   - :class:`covariance.EllipticEnvelop` is now deprecated - Please use :class:`covariance.EllipticEnvelope`
     instead.

   - ``NeighborsClassifier`` and ``NeighborsRegressor`` are gone in the module
     :ref:`neighbors`. Use the classes :class:`KNeighborsClassifier`,
     :class:`RadiusNeighborsClassifier`, :class:`KNeighborsRegressor`
     and/or :class:`RadiusNeighborsRegressor` instead.

   - Sparse classes in the :ref:`sgd` module are now deprecated.

   - In :class:`mixture.GMM`, :class:`mixture.DPGMM` and :class:`mixture.VBGMM`,
     parameters must be passed to an object when initialising it and not through
     ``fit``. Now ``fit`` will only accept the data as an input parameter.

   - methods ``rvs`` and ``decode`` in :class:`GMM` module are now deprecated.
     ``sample`` and ``score`` or ``predict`` should be used instead.

   - attribute ``_scores`` and ``_pvalues`` in univariate feature selection
     objects are now deprecated.
     ``scores_`` or ``pvalues_`` should be used instead.

   - In :class:`LogisticRegression`, :class:`LinearSVC`, :class:`SVC` and
     :class:`NuSVC`, the ``class_weight`` parameter is now an initialization
     parameter, not a parameter to fit. This makes grid searches
     over this parameter possible.

   - LFW ``data`` is now always shape ``(n_samples, n_features)`` to be
     consistent with the Olivetti faces dataset. Use ``images`` and
     ``pairs`` attribute to access the natural images shapes instead.

   - In :class:`svm.LinearSVC`, the meaning of the ``multi_class`` parameter
     changed.  Options now are ``'ovr'`` and ``'crammer_singer'``, with
     ``'ovr'`` being the default.  This does not change the default behavior
     but hopefully is less confusing.

   - Class :class:`feature_selection.text.Vectorizer` is deprecated and
     replaced by :class:`feature_selection.text.TfidfVectorizer`.

   - The preprocessor / analyzer nested structure for text feature
     extraction has been removed. All those features are
     now directly passed as flat constructor arguments
     to :class:`feature_selection.text.TfidfVectorizer` and
     :class:`feature_selection.text.CountVectorizer`, in particular the
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

   - Class :class:`feature_selection.text.TfidfVectorizer` now derives directly
     from :class:`feature_selection.text.CountVectorizer` to make grid
     search trivial.

   - methods ``rvs`` in :class:`_BaseHMM` module are now deprecated.
     ``sample`` should be used instead.

   - Beam pruning option in :class:`_BaseHMM` module is removed since it is
     difficult to be Cythonized. If you are interested, you can look in the
     history codes by git.

   - The SVMlight format loader now supports files with both zero-based and
     one-based column indices, since both occur "in the wild".

   - Arguments in class :class:`ShuffleSplit` are now consistent with
     :class:`StratifiedShuffleSplit`. Arguments ``test_fraction`` and
     ``train_fraction`` are deprecated and renamed to ``test_size`` and
     ``train_size`` and can accept both ``float`` and ``int``.

   - Arguments in class :class:`Bootstrap` are now consistent with
     :class:`StratifiedShuffleSplit`. Arguments ``n_test`` and
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

0.10
====

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

   - Fixed a memory leak in in :ref:`svm` module by `Brian Holt`_ (issue #367).

   - Faster tests by `Fabian Pedregosa`_ and others.

   - Silhouette Coefficient cluster analysis evaluation metric added as
     :func:`sklearn.metrics.silhouette_score` by Robert Layton.

   - Fixed a bug in :ref:`k_means` in the handling of the ``n_init`` parameter:
     the clustering algorithm used to be run ``n_init`` times but the last
     solution was retained instead of the best solution by `Olivier Grisel`_.

   - Minor refactoring in :ref:`sgd` module; consolidated dense and sparse
     predict methods; Enhanced test time performance by converting model
     parameters to fortran-style arrays after fitting (only multi-class).

   - Adjusted Mutual Information metric added as
     :func:`sklearn.metrics.adjusted_mutual_info_score` by Robert Layton.

   - Models like SVC/SVR/LinearSVC/LogisticRegression from libsvm/liblinear
     now support scaling of C regularization parameter by the number of
     samples by `Alexandre Gramfort`_.

   - New :ref:`Ensemble Methods <ensemble>` module by `Gilles Louppe`_ and
     `Brian Holt`_. The module comes with the random forest algorithm and the
     extra-trees method, along with documentation and examples.

   - :ref:`outlier_detection`: outlier and novelty detection, by
     `Virgile Fritsch`_.

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
     (:func:`sklearn.datasets.fetch_20newsgroups_vectorized`) by
     `Mathieu Blondel`_.

   - :ref:`multiclass` by `Lars Buitinck`_.

   - Utilities for fast computation of mean and variance for sparse matrices
     by `Mathieu Blondel`_.

   - Make :func:`sklearn.preprocessing.scale` and
     :class:`sklearn.preprocessing.Scaler` work on sparse matrices by
     `Olivier Grisel`_

   - Feature importances using decision trees and/or forest of trees,
     by `Gilles Louppe`_.

   - Parallel implementation of forests of randomized trees by
     `Gilles Louppe`_.

   - :class:`sklearn.cross_validation.ShuffleSplit` can subsample the train
     sets as well as the test sets by `Olivier Grisel`_.

   - Errors in the build of the documentation fixed by `Andreas Müller`_.


API changes summary
-------------------

Here are the code migration instructions when upgrading from scikit-learn
version 0.9:

  - Some estimators that may overwrite their inputs to save memory previously
    had ``overwrite_`` parameters; these have been replaced with ``copy_``
    parameters with exactly the opposite meaning.

    This particularly affects some of the estimators in :mod:`linear_model`.
    The default behavior is still to copy everything passed in.

  - The SVMlight dataset loader :func:`sklearn.datasets.load_svmlight_file` no
    longer supports loading two files at once; use ``load_svmlight_files``
    instead. Also, the (unused) ``buffer_mb`` parameter is gone.

  - Sparse estimators in the :ref:`sgd` module use dense parameter vector
    ``coef_`` instead of ``sparse_coef_``. This significantly improves
    test time performance.

  - The :ref:`covariance` module now has a robust estimator of
    covariance, the Minimum Covariance Determinant estimator.

  - Cluster evaluation metrics in :mod:`metrics.cluster` have been refactored
    but the changes are backwards compatible. They have been moved to the
    :mod:`metrics.cluster.supervised`, along with
    :mod:`metrics.cluster.unsupervised` which contains the Silhouette
    Coefficient.

  - The ``permutation_test_score`` function now behaves the same way as
    ``cross_val_score`` (i.e. uses the mean score across the folds.)

  - Cross Validation generators now use integer indices (``indices=True``)
    by default instead of boolean masks. This make it more intuitive to
    use with sparse matrix data.

  - The functions used for sparse coding, ``sparse_encode`` and
    ``sparse_encode_parallel`` have been combined into
    :func:`sklearn.decomposition.sparse_encode`, and the shapes of the arrays
    have been transposed for consistency with the matrix factorization setting,
    as opposed to the regression setting.

  - Fixed an off-by-one error in the SVMlight/LibSVM file format handling;
    files generated using :func:`sklearn.datasets.dump_svmlight_file` should be
    re-generated. (They should continue to work, but accidentally had one
    extra column of zeros prepended.)

  - ``BaseDictionaryLearning`` class replaced by ``SparseCodingMixin``.

  - :func:`sklearn.utils.extmath.fast_svd` has been renamed
    :func:`sklearn.utils.extmath.randomized_svd` and the default
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
   *  24  `Virgile Fritsch`_
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

0.9
===

scikit-learn 0.9 was released on September 2011, three months after the 0.8
release and includes the new modules :ref:`manifold`, :ref:`dirichlet_process`
as well as several new algorithms and documentation improvements.

This release also includes the dictionary-learning work developed by
`Vlad Niculae`_ as part of the `Google Summer of Code
<http://code.google.com/soc/>`_ program.



.. |banner1| image:: ./auto_examples/manifold/images/thumb/plot_compare_methods.png
   :target: auto_examples/manifold/plot_compare_methods.html

.. |banner2| image:: ./auto_examples/linear_model/images/thumb/plot_omp.png
   :target: auto_examples/linear_model/plot_omp.html

.. |banner3| image:: ./auto_examples/decomposition/images/thumb/plot_kernel_pca.png
   :target: auto_examples/decomposition/plot_kernel_pca.html

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
     and Python version thanks to Jean Kossaifi.

   - :ref:`Loader for libsvm/svmlight format <libsvm_loader>` by
     `Mathieu Blondel`_ and `Lars Buitinck`_

   - Documentation improvements: thumbnails in
     :ref:`example gallery <examples-index>` by `Fabian Pedregosa`_.

   - Important bugfixes in :ref:`svm` module (segfaults, bad
     performance) by `Fabian Pedregosa`_.

   - Added :ref:`multinomial_naive_bayes` and :ref:`bernoulli_naive_bayes`
     by `Lars Buitinck`_

   - Text feature extraction optimizations by Lars Buitinck

   - Chi-Square feature selection
     (:func:`feature_selection.univariate_selection.chi2`) by `Lars Buitinck`_.

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

   - Implementation of :class:`linear_model.LassoLarsCV`
     (cross-validated Lasso solver using the Lars algorithm) and
     :class:`linear_model.LassoLarsIC` (BIC/AIC model
     selection in Lars) by `Gael Varoquaux`_
     and `Alexandre Gramfort`_

   - Scalability improvements to :func:`metrics.roc_curve` by Olivier Hervieu

   - Distance helper functions :func:`metrics.pairwise.pairwise_distances`
     and :func:`metrics.pairwise.pairwise_kernels` by Robert Layton

   - :class:`Mini-Batch K-Means <cluster.MiniBatchKMeans>` by Nelle Varoquaux and Peter Prettenhofer.

   - :ref:`mldata` utilities by Pietro Berkes.

   - :ref:`olivetti_faces` by `David Warde-Farley`_.


API changes summary
-------------------

Here are the code migration instructions when upgrading from scikit-learn
version 0.8:

  - The ``scikits.learn`` package was renamed ``sklearn``. There is
    still a ``scikits.learn`` package alias for backward compatibility.

    Third-party projects with a dependency on scikit-learn 0.9+ should
    upgrade their codebase. For instance under Linux / MacOSX just run
    (make a backup first!)::

      find -name "*.py" | xargs sed -i 's/\bscikits.learn\b/sklearn/g'

  - Estimators no longer accept model parameters as ``fit`` arguments:
    instead all parameters must be only be passed as constructor
    arguments or using the now public ``set_params`` method inherited
    from :class:`base.BaseEstimator`.

    Some estimators can still accept keyword arguments on the ``fit``
    but this is restricted to data-dependent values (e.g. a Gram matrix
    or an affinity matrix that are precomputed from the ``X`` data matrix.

  - The ``cross_val`` package has been renamed to ``cross_validation``
    although there is also a ``cross_val`` package alias in place for
    backward compatibility.

    Third-party projects with a dependency on scikit-learn 0.9+ should
    upgrade their codebase. For instance under Linux / MacOSX just run
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
   - 32  Jean Kossaifi
   - 30  Conrad Lee
   - 22  Pietro Berkes
   - 18  andy
   - 17  David Warde-Farley
   - 12  Brian Holt
   - 11  Robert
   - 8  Amit Aides
   - 8  `Virgile Fritsch`_
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

0.8
===

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

  - :ref:`labeled_faces_in_the_wild` by `Olivier Grisel`_.

  - New :ref:`cross_decomposition` module by `Edouard Duchesnay`_.

  - :ref:`NMF` module `Vlad Niculae`_

  - Implementation of the :ref:`oracle_approximating_shrinkage` algorithm by
    `Virgile Fritsch`_ in the :ref:`covariance` module.


Some other modules benefited from significant improvements or cleanups.


  - Initial support for Python 3: builds and imports cleanly,
    some modules are usable while others have failing tests by `Fabian Pedregosa`_.

  - :class:`decomposition.PCA` is now usable from the Pipeline object by `Olivier Grisel`_.

  - Guide :ref:`performance-howto` by `Olivier Grisel`_.

  - Fixes for memory leaks in libsvm bindings, 64-bit safer BallTree by Lars Buitinck.

  - bug and style fixing in :ref:`k_means` algorithm by Jan Schlüter.

  - Add attribute converged to Gaussian Mixture Models by Vincent Schut.

  - Implemented ``transform``, ``predict_log_proba`` in :class:`lda.LDA`
    By `Mathieu Blondel`_.

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
   - 11  `Virgile Fritsch`_
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

0.7
===

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
    algorithm for finding barycenter weights. Also added some
    developer documentation for this module, see
    `notes_neighbors
    <https://github.com/scikit-learn/scikit-learn/wiki/Neighbors-working-notes>`_ for more information [`Fabian Pedregosa`_].

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

0.6
===

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
    sklearn.svm.sparse and sklearn.linear_model.sparse.

  - Lots of cool new examples and a new section that uses real-world
    datasets was created. These include:
    :ref:`example_applications_face_recognition.py`,
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

People that made this release possible preceded by number of commits:

   * 207  `Olivier Grisel`_

   * 167 `Fabian Pedregosa`_

   * 97 `Peter Prettenhofer`_

   * 68 `Alexandre Gramfort`_

   * 59  `Mathieu Blondel`_

   * 55  `Gael Varoquaux`_

   * 33  Vincent Dubourg

   * 21  `Ron Weiss <http://www.ee.columbia.edu/~ronw/>`_

   * 9  Bertrand Thirion

   * 3  `Alexandre Passos`_

   * 3  Anne-Laure Fouque

   * 2  Ronan Amicel

   * 1 `Christian Osendorfer`_



.. _changes_0_5:


0.5
===

Changelog
---------

New classes
-----------

    - Support for sparse matrices in some classifiers of modules
      ``svm`` and ``linear_model`` (see :class:`svm.sparse.SVC`,
      :class:`svm.sparse.SVR`, :class:`svm.sparse.LinearSVC`,
      :class:`linear_model.sparse.Lasso`, :class:`linear_model.sparse.ElasticNet`)

    - New :class:`pipeline.Pipeline` object to compose different estimators.

    - Recursive Feature Elimination routines in module
      :ref:`feature_selection`.

    - Addition of various classes capable of cross validation in the
      linear_model module (:class:`linear_model.LassoCV`, :class:`linear_model.ElasticNetCV`,
      etc.).

    - New, more efficient LARS algorithm implementation. The Lasso
      variant of the algorithm is also implemented. See
      :class:`linear_model.lars_path`, :class:`linear_model.Lars` and
      :class:`linear_model.LassoLars`.

    - New Hidden Markov Models module (see classes
      :class:`hmm.GaussianHMM`, :class:`hmm.MultinomialHMM`,
      :class:`hmm.GMMHMM`)

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
      ``example_mlcomp_sparse_document_classification.py`` (since removed) and
      :ref:`example_text_document_classification_20newsgroups.py`

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

.. _Olivier Grisel: http://twitter.com/ogrisel

.. _Gael Varoquaux: http://gael-varoquaux.info

.. _Alexandre Gramfort: http://alexandre.gramfort.net

.. _Fabian Pedregosa: http://fseoane.net/blog/

.. _Mathieu Blondel: http://www.mblondel.org

.. _James Bergstra: http://www-etud.iro.umontreal.ca/~bergstrj/

.. _liblinear: http://www.csie.ntu.edu.tw/~cjlin/liblinear/

.. _Yaroslav Halchenko: http://www.onerussian.com/

.. _Vlad Niculae: http://vene.ro

.. _Edouard Duchesnay: http://www.lnao.fr/spip.php?rubrique30

.. _Peter Prettenhofer: http://sites.google.com/site/peterprettenhofer/

.. _Alexandre Passos: <http://atpassos.posterous.com>

.. _Nicolas Pinto: http://pinto.scripts.mit.edu/

.. _Virgile Fritsch: http://parietal.saclay.inria.fr/Members/virgile-fritsch

.. _Bertrand Thirion: http://parietal.saclay.inria.fr/Members/bertrand-thirion

.. _Andreas Müller: http://peekaboo-vision.blogspot.com

.. _Matthieu Perrot: http://www.lnao.fr/spip.php?rubrique19

.. _Jake Vanderplas: http://www.astro.washington.edu/users/vanderplas/

.. _Gilles Louppe: http://www.montefiore.ulg.ac.be/~glouppe/

.. _INRIA: http://inria.fr

.. _Parietal Team: http://parietal.saclay.inria.fr/

.. _Lars Buitinck: https://github.com/larsmans

.. _David Warde-Farley: http://www-etud.iro.umontreal.ca/~wardefar/

.. _Brian Holt: http://info.ee.surrey.ac.uk/Personal/B.Holt/

.. _Satrajit Ghosh: http://www.mit.edu/~satra/

.. _Robert Layton: http://www.twitter.com/robertlayton

.. _Scott White: http://twitter.com/scottblanc

.. _Jaques Grobler: https://github.com/jaquesgrobler/scikit-learn/wiki/Jaques-Grobler

.. _David Marek: http://www.davidmarek.cz/

.. _@kernc: http://github.com/kernc

.. _Christian Osendorfer: http://osdf.github.com

.. _Noel Dawe: http://noel.dawe.me

.. _Arnaud Joly: http://www.ajoly.org

.. _Rob Zinkov: http://zinkov.com

.. _Martin Luessi: https://github.com/mluessi

.. _Joel Nothman: http://joelnothman.com

.. _Norbert Crombach: https://github.com/norbert

.. _Eustache Diemert: https://github.com/oddskool

.. _Justin Vincent: https://github.com/justinvf

.. _Denis Engemann: https://github.com/dengemann

.. _Nicolas Trésegnie : http://nicolastr.com/

.. _Kemal Eren: http://www.kemaleren.com

.. _Yann Dauphin: http://ynd.github.io/

.. _Nelle Varoquaux: https://github.com/nellev

.. _Subhodeep Moitra: https://github.com/smoitra87

.. _Yannick Schwartz: https://team.inria.fr/parietal/schwarty/

.. _Mikhail Korobov: http://kmike.ru/pages/about/

.. _Kyle Kastner: http://kastnerkyle.github.io

.. _@FedericoV: https://github.com/FedericoV/

.. _Daniel Nouri: http://danielnouri.org

.. _Johannes Schönberger: https://github.com/ahojnnes

.. _Manoj Kumar: https://manojbits.wordpress.com

.. _Andrew Tulloch: http://tullo.ch

.. _Maheshakya Wijewardena: https://github.com/maheshakya

.. _Danny Sullivan: https://github.com/dsullivan7

.. _Michael Eickenberg: https://github.com/eickenberg

.. _Jeffrey Blackburne: https://github.com/jblackburne

.. _Hamzeh Alsalhi: https://github.com/hamsal

.. _Ronald Phlypo: https://github.com/rphlypo

.. _Laurent Direr: https://github.com/ldirer

.. _Nikolay Mayorov: https://github.com/nmayorov

.. _Jatin Shah: http://jatinshah.org/

.. _Dougal Sutherland: https://github.com/dougalsutherland

.. _Michal Romaniuk: https://github.com/romaniukm
