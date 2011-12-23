=========
Reference
=========

This is the class and function reference of scikit-learn. Please refer to
the :ref:`full user guide <user_guide>` for further details, as the class and
function raw specifications may not be enough to give full guidelines on their
uses.

.. contents:: List of modules
   :local:


.. _cluster_ref:

:mod:`sklearn.cluster`: Clustering
==================================

.. automodule:: sklearn.cluster
   :no-members:
   :no-inherited-members:

**User guide:** See the :ref:`clustering` section for further details.

.. currentmodule:: sklearn

.. autosummary::
   :toctree: generated/
   :template: class.rst

   cluster.AffinityPropagation
   cluster.DBSCAN
   cluster.KMeans
   cluster.MiniBatchKMeans
   cluster.MeanShift
   cluster.SpectralClustering
   cluster.Ward


.. _covariance_ref:

:mod:`sklearn.covariance`: Covariance Estimators
================================================

.. automodule:: sklearn.covariance
   :no-members:
   :no-inherited-members:

**User guide:** See the :ref:`covariance` section for further details.

.. currentmodule:: sklearn

.. autosummary::
   :toctree: generated/
   :template: class.rst

   covariance.EmpiricalCovariance
   covariance.ShrunkCovariance
   covariance.LedoitWolf
   covariance.OAS
   covariance.GraphLasso
   covariance.GraphLassoCV
   covariance.MinCovDet

.. autosummary::
   :toctree: generated/
   :template: function.rst

   covariance.empirical_covariance
   covariance.ledoit_wolf
   covariance.shrunk_covariance
   covariance.oas
   covariance.graph_lasso


.. _cross_validation_ref:

:mod:`sklearn.cross_validation`: Cross Validation
=================================================

.. automodule:: sklearn.cross_validation
   :no-members:
   :no-inherited-members:

**User guide:** See the :ref:`cross_validation` section for further details.

.. currentmodule:: sklearn

.. autosummary::
   :toctree: generated/
   :template: class.rst

   cross_validation.LeaveOneOut
   cross_validation.LeavePOut
   cross_validation.KFold
   cross_validation.StratifiedKFold
   cross_validation.LeaveOneLabelOut
   cross_validation.LeavePLabelOut
   cross_validation.Bootstrap
   cross_validation.ShuffleSplit


.. _datasets_ref:

:mod:`sklearn.datasets`: Datasets
=================================

.. automodule:: sklearn.datasets
   :no-members:
   :no-inherited-members:

**User guide:** See the :ref:`datasets` section for further details.

Loaders
-------

.. currentmodule:: sklearn

.. autosummary::
   :toctree: generated/
   :template: function.rst

   datasets.load_boston
   datasets.load_files
   datasets.load_diabetes
   datasets.load_digits
   datasets.load_iris
   datasets.load_linnerud
   datasets.load_lfw_pairs
   datasets.fetch_lfw_pairs
   datasets.load_lfw_people
   datasets.fetch_lfw_people
   datasets.load_20newsgroups
   datasets.fetch_20newsgroups
   datasets.fetch_20newsgroups_vectorized
   datasets.fetch_olivetti_faces

Samples generator
-----------------

.. currentmodule:: sklearn

.. autosummary::
   :toctree: generated/
   :template: function.rst

   datasets.make_classification
   datasets.make_multilabel_classification
   datasets.make_regression
   datasets.make_blobs
   datasets.make_friedman1
   datasets.make_friedman2
   datasets.make_friedman3
   datasets.make_low_rank_matrix
   datasets.make_sparse_coded_signal
   datasets.make_sparse_uncorrelated
   datasets.make_spd_matrix
   datasets.make_sparse_spd_matrix
   datasets.make_swiss_roll
   datasets.make_s_curve


.. _decomposition_ref:

:mod:`sklearn.decomposition`: Matrix Decomposition
==================================================

.. automodule:: sklearn.decomposition
   :no-members:
   :no-inherited-members:

**User guide:** See the :ref:`decompositions` section for further details.

.. currentmodule:: sklearn

.. autosummary::
   :toctree: generated/
   :template: class.rst

   decomposition.PCA
   decomposition.ProbabilisticPCA
   decomposition.ProjectedGradientNMF
   decomposition.RandomizedPCA
   decomposition.KernelPCA
   decomposition.FastICA
   decomposition.NMF
   decomposition.SparsePCA
   decomposition.MiniBatchSparsePCA
   decomposition.SparseCoder
   decomposition.DictionaryLearning
   decomposition.MiniBatchDictionaryLearning

.. autosummary::
   :toctree: generated/
   :template: function.rst

   decomposition.fastica
   decomposition.dict_learning
   decomposition.dict_learning_online
   decomposition.sparse_encode


.. _ensemble_ref:

:mod:`sklearn.ensemble`: Ensemble Methods
=========================================

.. automodule:: sklearn.ensemble
   :no-members:
   :no-inherited-members:

**User guide:** See the :ref:`ensemble` section for further details.

.. currentmodule:: sklearn

.. autosummary::
   :toctree: generated/
   :template: class.rst

   ensemble.RandomForestClassifier
   ensemble.RandomForestRegressor
   ensemble.ExtraTreesClassifier
   ensemble.ExtraTreesRegressor

.. autosummary::
   :toctree: generated/
   :template: function.rst


.. _feature_extraction_ref:

:mod:`sklearn.feature_extraction`: Feature Extraction
=====================================================

.. automodule:: sklearn.feature_extraction
   :no-members:
   :no-inherited-members:

**User guide:** See the :ref:`feature_extraction` section for further details.

.. currentmodule:: sklearn

From images
-----------

.. automodule:: sklearn.feature_extraction.image
   :no-members:
   :no-inherited-members:

.. currentmodule:: sklearn

.. autosummary::
   :toctree: generated/
   :template: function.rst

   feature_extraction.image.img_to_graph
   feature_extraction.image.grid_to_graph
   feature_extraction.image.extract_patches_2d
   feature_extraction.image.reconstruct_from_patches_2d

   :template: class.rst

   feature_extraction.image.PatchExtractor

From text
---------

.. automodule:: sklearn.feature_extraction.text
   :no-members:
   :no-inherited-members:

.. currentmodule:: sklearn

.. autosummary::
   :toctree: generated/
   :template: class.rst

   feature_extraction.text.RomanPreprocessor
   feature_extraction.text.WordNGramAnalyzer
   feature_extraction.text.CharNGramAnalyzer
   feature_extraction.text.CountVectorizer
   feature_extraction.text.TfidfTransformer
   feature_extraction.text.Vectorizer


.. _feature_selection_ref:

:mod:`sklearn.feature_selection`: Feature Selection
===================================================

.. automodule:: sklearn.feature_selection
   :no-members:
   :no-inherited-members:

**User guide:** See the :ref:`feature_selection` section for further details.

.. currentmodule:: sklearn

.. autosummary::
   :toctree: generated/
   :template: class.rst

   feature_selection.SelectPercentile
   feature_selection.SelectKBest
   feature_selection.SelectFpr
   feature_selection.SelectFdr
   feature_selection.SelectFwe
   feature_selection.RFE
   feature_selection.RFECV

.. autosummary::
   :toctree: generated/
   :template: function.rst

   feature_selection.chi2
   feature_selection.f_classif
   feature_selection.f_regression


.. _gaussian_process_ref:

:mod:`sklearn.gaussian_process`: Gaussian Processes
===================================================

.. automodule:: sklearn.gaussian_process
   :no-members:
   :no-inherited-members:

**User guide:** See the :ref:`gaussian_process` section for further details.

.. currentmodule:: sklearn

.. autosummary::
  :toctree: generated/
  :template: class.rst

  gaussian_process.GaussianProcess

.. autosummary::
   :toctree: generated
   :template: function.rst

   gaussian_process.correlation_models.absolute_exponential
   gaussian_process.correlation_models.squared_exponential
   gaussian_process.correlation_models.generalized_exponential
   gaussian_process.correlation_models.pure_nugget
   gaussian_process.correlation_models.cubic
   gaussian_process.correlation_models.linear
   gaussian_process.regression_models.constant
   gaussian_process.regression_models.linear
   gaussian_process.regression_models.quadratic


.. _grid_search_ref:

:mod:`sklearn.grid_search`: Grid Search
=======================================

.. automodule:: sklearn.grid_search
   :no-members:
   :no-inherited-members:

**User guide:** See the :ref:`grid_search` section for further details.

.. currentmodule:: sklearn

.. autosummary::
   :toctree: generated/
   :template: class.rst

   grid_search.GridSearchCV
   grid_search.IterGrid


.. _hmm_ref:

:mod:`sklearn.hmm`: Hidden Markov Models
========================================

.. automodule:: sklearn.hmm
   :no-members:
   :no-inherited-members:

**User guide:** See the :ref:`hmm` section for further details.

.. currentmodule:: sklearn

.. autosummary::
   :toctree: generated/
   :template: class.rst

   hmm.GaussianHMM
   hmm.MultinomialHMM
   hmm.GMMHMM


.. _kernel_approximation_ref:

:mod:`sklearn.kernel_approximation` Kernel Approximation
========================================================

.. automodule:: sklearn.kernel_approximation
   :no-members:
   :no-inherited-members:

**User guide:** See the :ref:`kernel_approximation` section for further details.

.. currentmodule:: sklearn

.. autosummary::
   :toctree: generated/
   :template: class.rst

   kernel_approximation.RBFSampler
   kernel_approximation.AdditiveChi2Sampler
   kernel_approximation.SkewedChi2Sampler


.. _lda_ref:

:mod:`sklearn.lda`: Linear Discriminant Analysis
================================================

.. automodule:: sklearn.lda
   :no-members:
   :no-inherited-members:

.. currentmodule:: sklearn

.. autosummary::
   :toctree: generated
   :template: class.rst

   lda.LDA


.. _linear_model_ref:

:mod:`sklearn.linear_model`: Generalized Linear Models
======================================================

.. automodule:: sklearn.linear_model
   :no-members:
   :no-inherited-members:

**User guide:** See the :ref:`linear_model` section for further details.

For dense data
--------------

.. currentmodule:: sklearn

.. autosummary::
   :toctree: generated/
   :template: class.rst

   linear_model.LinearRegression
   linear_model.Ridge
   linear_model.RidgeCV
   linear_model.Lasso
   linear_model.LassoCV
   linear_model.ElasticNet
   linear_model.ElasticNetCV
   linear_model.Lars
   linear_model.LassoLars
   linear_model.LarsCV
   linear_model.LassoLarsCV
   linear_model.LassoLarsIC
   linear_model.LogisticRegression
   linear_model.OrthogonalMatchingPursuit
   linear_model.SGDClassifier
   linear_model.SGDRegressor
   linear_model.BayesianRidge
   linear_model.ARDRegression

.. autosummary::
   :toctree: generated/
   :template: function.rst

   linear_model.lasso_path
   linear_model.lars_path
   linear_model.orthogonal_mp
   linear_model.orthogonal_mp_gram

For sparse data
---------------

.. automodule:: sklearn.linear_model.sparse
   :no-members:
   :no-inherited-members:

**User guide:** See the :ref:`linear_model` section for further details.

.. currentmodule:: sklearn

.. autosummary::
   :toctree: generated/
   :template: class.rst

   linear_model.sparse.Lasso
   linear_model.sparse.ElasticNet
   linear_model.sparse.SGDClassifier
   linear_model.sparse.SGDRegressor
   linear_model.sparse.LogisticRegression


.. _manifold_ref:

:mod:`sklearn.manifold`: Manifold Learning
==========================================

.. automodule:: sklearn.manifold
   :no-members:
   :no-inherited-members:

**User guide:** See the :ref:`manifold` section for further details.

.. currentmodule:: sklearn

.. autosummary::
    :toctree: generated
    :template: class.rst

    manifold.LocallyLinearEmbedding
    manifold.Isomap

.. autosummary::
    :toctree: generated
    :template: function.rst

    manifold.locally_linear_embedding


.. _metrics_ref:

:mod:`sklearn.metrics`: Metrics
===============================

.. automodule:: sklearn.metrics
   :no-members:
   :no-inherited-members:

.. currentmodule:: sklearn

Classification metrics
----------------------

.. autosummary::
   :toctree: generated/
   :template: function.rst

   metrics.confusion_matrix
   metrics.roc_curve
   metrics.auc
   metrics.precision_score
   metrics.recall_score
   metrics.fbeta_score
   metrics.f1_score
   metrics.precision_recall_fscore_support
   metrics.classification_report
   metrics.precision_recall_curve
   metrics.zero_one_score
   metrics.zero_one
   metrics.hinge_loss

Regression metrics
------------------

.. autosummary::
   :toctree: generated/
   :template: function.rst

   metrics.r2_score
   metrics.mean_square_error

Clustering metrics
------------------

See the :ref:`clustering` section of the user guide for further details.

.. automodule:: sklearn.metrics.cluster
   :no-members:
   :no-inherited-members:

.. currentmodule:: sklearn

.. autosummary::
   :toctree: generated/
   :template: function.rst

   metrics.adjusted_rand_score
   metrics.adjusted_mutual_info_score
   metrics.homogeneity_completeness_v_measure
   metrics.homogeneity_score
   metrics.completeness_score
   metrics.v_measure_score
   metrics.silhouette_score

Pairwise metrics
----------------

.. automodule:: sklearn.metrics.pairwise
   :no-members:
   :no-inherited-members:

.. currentmodule:: sklearn

.. autosummary::
   :toctree: generated/
   :template: function.rst

   metrics.pairwise.euclidean_distances
   metrics.pairwise.manhattan_distances
   metrics.pairwise.linear_kernel
   metrics.pairwise.polynomial_kernel
   metrics.pairwise.rbf_kernel
   metrics.pairwise.distance_metrics
   metrics.pairwise.pairwise_distances
   metrics.pairwise.kernel_metrics
   metrics.pairwise.pairwise_kernels


.. _mixture_ref:

:mod:`sklearn.mixture`: Gaussian Mixture Models
===============================================

.. automodule:: sklearn.mixture
   :no-members:
   :no-inherited-members:

**User guide:** See the :ref:`mixture` section for further details.

.. currentmodule:: sklearn

.. autosummary::
   :toctree: generated/
   :template: class.rst

   mixture.GMM
   mixture.DPGMM
   mixture.VBGMM


.. _multiclass_ref:

:mod:`sklearn.multiclass`: Multiclass and multilabel classification
===================================================================

.. automodule:: sklearn.multiclass
   :no-members:
   :no-inherited-members:

**User guide:** See the :ref:`multiclass` section for further details.

.. currentmodule:: sklearn

.. autosummary::
    :toctree: generated
    :template: class.rst

    multiclass.OneVsRestClassifier
    multiclass.OneVsOneClassifier
    multiclass.OutputCodeClassifier

.. autosummary::
    :toctree: generated
    :template: function.rst

    multiclass.fit_ovr
    multiclass.predict_ovr
    multiclass.fit_ovo
    multiclass.predict_ovo
    multiclass.fit_ecoc
    multiclass.predict_ecoc


.. _naive_bayes_ref:

:mod:`sklearn.naive_bayes`: Naive Bayes
=======================================

.. automodule:: sklearn.naive_bayes
   :no-members:
   :no-inherited-members:

**User guide:** See the :ref:`naive_bayes` section for further details.

.. currentmodule:: sklearn

.. autosummary::
   :toctree: generated/
   :template: class.rst

   naive_bayes.GaussianNB
   naive_bayes.MultinomialNB
   naive_bayes.BernoulliNB


.. _neighbors_ref:

:mod:`sklearn.neighbors`: Nearest Neighbors
===========================================

.. automodule:: sklearn.neighbors
   :no-members:
   :no-inherited-members:

**User guide:** See the :ref:`neighbors` section for further details.

.. currentmodule:: sklearn

.. autosummary::
   :toctree: generated/
   :template: class.rst

   neighbors.NearestNeighbors
   neighbors.KNeighborsClassifier
   neighbors.RadiusNeighborsClassifier
   neighbors.NeighborsClassifier
   neighbors.KNeighborsRegressor
   neighbors.RadiusNeighborsRegressor
   neighbors.NeighborsRegressor
   neighbors.BallTree

.. autosummary::
   :toctree: generated/
   :template: function.rst

   neighbors.kneighbors_graph
   neighbors.radius_neighbors_graph


.. _pls_ref:

:mod:`sklearn.pls`: Partial Least Squares
=========================================

.. automodule:: sklearn.pls
   :no-members:
   :no-inherited-members:

**User guide:** See the :ref:`pls` section for further details.

.. currentmodule:: sklearn

.. autosummary::
   :toctree: generated/
   :template: class.rst

   pls.PLSRegression
   pls.PLSCanonical
   pls.CCA
   pls.PLSSVD


.. _pipeline_ref:

:mod:`sklearn.pipeline`: Pipeline
=================================

.. automodule:: sklearn.pipeline
   :no-members:
   :no-inherited-members:

.. currentmodule:: sklearn

.. autosummary::
   :toctree: generated/
   :template: class.rst

   pipeline.Pipeline


.. _preprocessing_ref:

:mod:`sklearn.preprocessing`: Preprocessing and Normalization
=============================================================

.. automodule:: sklearn.preprocessing
   :no-members:
   :no-inherited-members:

**User guide:** See the :ref:`preprocessing` section for further details.

.. currentmodule:: sklearn

.. autosummary::
   :toctree: generated/
   :template: class.rst

   preprocessing.Scaler
   preprocessing.Normalizer
   preprocessing.Binarizer
   preprocessing.LabelBinarizer
   preprocessing.KernelCenterer

.. autosummary::
   :toctree: generated/
   :template: function.rst

   preprocessing.scale
   preprocessing.normalize
   preprocessing.binarize


.. _svm_ref:

:mod:`sklearn.svm`: Support Vector Machines
===========================================

.. automodule:: sklearn.svm
   :no-members:
   :no-inherited-members:

**User guide:** See the :ref:`svm` section for further details.

For dense data
--------------

.. currentmodule:: sklearn

.. autosummary::
   :toctree: generated/
   :template: class.rst

   svm.SVC
   svm.LinearSVC
   svm.NuSVC
   svm.SVR
   svm.NuSVR
   svm.OneClassSVM

.. autosummary::
   :toctree: generated/
   :template: function.rst

   svm.l1_min_c

For sparse data
---------------

.. automodule:: sklearn.svm.sparse
   :no-members:
   :no-inherited-members:

.. currentmodule:: sklearn

.. autosummary::
   :toctree: generated/
   :template: class.rst

   svm.sparse.SVC
   svm.sparse.NuSVC
   svm.sparse.SVR
   svm.sparse.NuSVR
   svm.sparse.OneClassSVM
   svm.sparse.LinearSVC

Low-level methods
-----------------

.. autosummary::
   :toctree: generated
   :template: function.rst

   svm.libsvm.fit
   svm.libsvm.decision_function
   svm.libsvm.predict
   svm.libsvm.predict_proba
   svm.libsvm.cross_validation


.. _tree_ref:

:mod:`sklearn.tree`: Decision Trees
===================================

.. automodule:: sklearn.tree
   :no-members:
   :no-inherited-members:

**User guide:** See the :ref:`tree` section for further details.

.. currentmodule:: sklearn

.. autosummary::
   :toctree: generated/
   :template: class.rst

   tree.DecisionTreeClassifier
   tree.DecisionTreeRegressor
   tree.ExtraTreeClassifier
   tree.ExtraTreeRegressor

.. autosummary::
   :toctree: generated/
   :template: function.rst

   tree.export_graphviz


.. _utils_ref:

:mod:`sklearn.utils`: Utilities
===============================

.. automodule:: sklearn.utils
   :no-members:
   :no-inherited-members:

**Developer guide:** See the :ref:`developers-utils` page for further details.

.. currentmodule:: sklearn

.. autosummary::
   :toctree: generated/
   :template: function.rst

   utils.check_random_state
   utils.resample
   utils.shuffle
