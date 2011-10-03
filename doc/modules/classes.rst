===============
Class reference
===============


.. TODO: put some order here. Alphabetical ?


Support Vector Machines
=======================


.. automodule:: sklearn.svm
   :no-members:
   :no-inherited-members:

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


Generalized Linear Models
=========================

.. automodule:: sklearn.linear_model
   :no-members:
   :no-inherited-members:

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

.. currentmodule:: sklearn

.. autosummary::
   :toctree: generated/
   :template: class.rst

   linear_model.sparse.Lasso
   linear_model.sparse.ElasticNet
   linear_model.sparse.SGDClassifier
   linear_model.sparse.SGDRegressor
   linear_model.sparse.LogisticRegression


Naive Bayes
===========

.. automodule:: sklearn.naive_bayes
   :no-members:
   :no-inherited-members:

.. currentmodule:: sklearn

.. autosummary::
   :toctree: generated/
   :template: class.rst

   naive_bayes.GaussianNB
   naive_bayes.MultinomialNB
   naive_bayes.BernoulliNB


Nearest Neighbors
=================

.. automodule:: sklearn.neighbors
   :no-members:
   :no-inherited-members:

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


Gaussian Mixture Models
=======================

.. automodule:: sklearn.mixture
   :no-members:
   :no-inherited-members:

.. currentmodule:: sklearn

.. autosummary::
   :toctree: generated/
   :template: class.rst

   mixture.GMM
   mixture.DPGMM
   mixture.VBGMM


Hidden Markov Models
====================

.. automodule:: sklearn.hmm
   :no-members:
   :no-inherited-members:

.. currentmodule:: sklearn

.. autosummary::
   :toctree: generated/
   :template: class.rst

   hmm.GaussianHMM
   hmm.MultinomialHMM
   hmm.GMMHMM


Clustering
==========

.. automodule:: sklearn.cluster
   :no-members:
   :no-inherited-members:

.. currentmodule:: sklearn

.. autosummary::
   :toctree: generated/
   :template: class.rst

   cluster.KMeans
   cluster.MiniBatchKMeans
   cluster.MeanShift
   cluster.SpectralClustering
   cluster.AffinityPropagation
   cluster.Ward


Metrics
=======


Classification metrics
----------------------

.. automodule:: sklearn.metrics
   :no-members:
   :no-inherited-members:

.. currentmodule:: sklearn

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

.. automodule:: sklearn.metrics
   :no-members:
   :no-inherited-members:

.. currentmodule:: sklearn

.. autosummary::
   :toctree: generated/
   :template: function.rst

   metrics.r2_score
   metrics.mean_square_error


Clustering metrics
------------------

.. automodule:: sklearn.metrics.cluster
   :no-members:
   :no-inherited-members:

.. currentmodule:: sklearn

.. autosummary::
   :toctree: generated/
   :template: function.rst

   metrics.adjusted_rand_score
   metrics.homogeneity_completeness_v_measure
   metrics.homogeneity_score
   metrics.completeness_score
   metrics.v_measure_score


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
   metrics.pairwise.l1_distances
   metrics.pairwise.linear_kernel
   metrics.pairwise.polynomial_kernel
   metrics.pairwise.rbf_kernel
   metrics.pairwise.distance_metrics
   metrics.pairwise.pairwise_distances
   metrics.pairwise.kernel_metrics
   metrics.pairwise.pairwise_kernels

Covariance Estimators
=====================

.. automodule:: sklearn.covariance
   :no-members:
   :no-inherited-members:

.. currentmodule:: sklearn

.. autosummary::
   :toctree: generated/
   :template: class.rst

   covariance.EmpiricalCovariance
   covariance.ShrunkCovariance
   covariance.LedoitWolf
   covariance.OAS

.. autosummary::
   :toctree: generated/
   :template: function.rst

   covariance.empirical_covariance
   covariance.ledoit_wolf
   covariance.shrunk_covariance
   covariance.oas


Signal Decomposition
====================

.. automodule:: sklearn.decomposition
   :no-members:
   :no-inherited-members:

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
   decomposition.DictionaryLearning
   decomposition.MiniBatchDictionaryLearning

.. autosummary::
   :toctree: generated/
   :template: function.rst
 
   decomposition.fastica
   decomposition.dict_learning
   decomposition.dict_learning_online
   decomposition.sparse_encode
   decomposition.sparse_encode_parallel


Linear Discriminant Analysis
============================

.. autosummary::
   :toctree: generated
   :template: class.rst

   lda.LDA


Partial Least Squares
=====================

.. automodule:: sklearn.pls
   :no-members:
   :no-inherited-members:

.. currentmodule:: sklearn

.. autosummary::
   :toctree: generated/
   :template: class.rst

   pls.PLSRegression
   pls.PLSCanonical
   pls.CCA
   pls.PLSSVD


Cross Validation
================

.. automodule:: sklearn.cross_validation
   :no-members:
   :no-inherited-members:

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


Grid Search
===========

.. automodule:: sklearn.grid_search
   :no-members:
   :no-inherited-members:

.. currentmodule:: sklearn

.. autosummary::
   :toctree: generated/
   :template: class.rst

   grid_search.GridSearchCV


.. _feature_selection_ref:


Feature Selection
=================

.. automodule:: sklearn.feature_selection
   :no-members:
   :no-inherited-members:

.. currentmodule:: sklearn

.. autosummary::
   :toctree: generated/
   :template: class.rst

   feature_selection.rfe.RFE
   feature_selection.rfe.RFECV


.. autosummary::
   :toctree: generated/
   :template: function.rst

   feature_selection.univariate_selection.chi2


.. _feature_extraction:


Feature Extraction
==================

.. automodule:: sklearn.feature_extraction
   :no-members:
   :no-inherited-members:

.. currentmodule:: sklearn

From images
-------------

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


Preprocessing and normalization
===============================

.. automodule:: sklearn.preprocessing
   :no-members:
   :no-inherited-members:

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


Manifold learning
=================

.. autosummary::
    :toctree: generated
    :template: class.rst

    manifold.LocallyLinearEmbedding
    manifold.Isomap


.. autosummary::
    :toctree: generated
    :template: function.rst

    manifold.locally_linear_embedding


Datasets
========

Loaders
-------

.. automodule:: sklearn.datasets
   :no-members:
   :no-inherited-members:

.. currentmodule:: sklearn

.. autosummary::
   :toctree: generated/
   :template: function.rst

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

Samples generator
-----------------

.. automodule:: sklearn.datasets
   :no-members:
   :no-inherited-members:

.. currentmodule:: sklearn

.. autosummary::
   :toctree: generated/
   :template: function.rst

   datasets.make_classification
   datasets.make_regression
   datasets.make_blobs
   datasets.make_friedman1
   datasets.make_friedman2
   datasets.make_friedman3
   datasets.make_low_rank_matrix
   datasets.make_sparse_coded_signal
   datasets.make_sparse_uncorrelated
   datasets.make_spd_matrix
   datasets.make_swiss_roll
   datasets.make_s_curve


Pipeline
========

.. automodule:: sklearn.pipeline
   :no-members:
   :no-inherited-members:

.. currentmodule:: sklearn

.. autosummary::
   :toctree: generated/
   :template: class.rst

   pipeline.Pipeline


Utilities
=========

.. automodule:: sklearn.utils
   :no-members:
   :no-inherited-members:

.. currentmodule:: sklearn

.. autosummary::
   :toctree: generated/
   :template: function.rst

   utils.check_random_state
   utils.resample
   utils.shuffle
