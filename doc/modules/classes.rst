===============
Class reference
===============


.. TODO: put some order here. Alphabetical ?


Support Vector Machines
=======================

.. currentmodule:: scikits.learn

.. autosummary::
   :toctree: generated/
   :template: class.rst

   svm.SVC
   svm.LinearSVC
   svm.NuSVC
   svm.SVR
   svm.NuSVR
   svm.OneClassSVM

For sparse data
---------------

.. autosummary::
   :toctree: generated/
   :template: class.rst

   svm.sparse.SVC
   svm.sparse.NuSVC
   svm.sparse.SVR
   svm.sparse.NuSVR
   svm.sparse.OneClassSVM
   svm.sparse.LinearSVC


Generalized Linear Models
=========================

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
   linear_model.LARS
   linear_model.LassoLARS
   linear_model.LogisticRegression
   linear_model.SGDClassifier
   linear_model.SGDRegressor


.. autosummary::
   :toctree: generated/
   :template: function.rst

   linear_model.lasso_path
   linear_model.lars_path


For sparse data
---------------

.. autosummary::
   :toctree: generated/
   :template: class.rst

   linear_model.sparse.Lasso
   linear_model.sparse.ElasticNet
   linear_model.sparse.SGDClassifier
   linear_model.sparse.SGDRegressor


Bayesian Regression
===================

.. autosummary::
   :toctree: generated/
   :template: class.rst

   linear_model.BayesianRidge
   linear_model.ARDRegression


Naive Bayes
===========

.. autosummary::
   :toctree: generated/
   :template: class.rst

   naive_bayes.GNB


Nearest Neighbors
=================

.. autosummary::
   :toctree: generated/
   :template: class.rst

   neighbors.NeighborsClassifier
   neighbors.NeighborsRegressor
   ball_tree.BallTree

.. autosummary::
   :toctree: generated/
   :template: function.rst

   neighbors.kneighbors_graph


Gaussian Mixture Models
=======================

.. autosummary::
   :toctree: generated/
   :template: class.rst

   mixture.GMM


Hidden Markov Models
====================

.. autosummary::
   :toctree: generated/
   :template: class.rst

   hmm.GaussianHMM
   hmm.MultinomialHMM
   hmm.GMMHMM


Clustering
==========

.. autosummary::
   :toctree: generated/
   :template: class.rst

   cluster.KMeans
   cluster.MeanShift
   cluster.SpectralClustering
   cluster.AffinityPropagation
   cluster.Ward


Metrics
=======


.. autosummary::
   :toctree: generated/
   :template: function.rst

   metrics.euclidean_distances
   metrics.unique_labels
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
   metrics.r2_score
   metrics.zero_one_score
   metrics.zero_one
   metrics.mean_square_error


Covariance Estimators
=====================

.. autosummary::
   :toctree: generated/
   :template: class.rst

   covariance.Covariance
   covariance.ShrunkCovariance
   covariance.LedoitWolf

.. autosummary::
   :toctree: generated/
   :template: function.rst

   covariance.ledoit_wolf


Signal Decomposition
====================

.. autosummary::
   :toctree: generated/
   :template: class.rst

   pca.PCA
   pca.ProbabilisticPCA
   pca.RandomizedPCA
   fastica.FastICA
   nmf.NMF

.. autosummary::
   :toctree: generated/
   :template: function.rst

   fastica.fastica

Cross Validation
================

.. autosummary::
   :toctree: generated/
   :template: class.rst

   cross_val.LeaveOneOut
   cross_val.LeavePOut
   cross_val.KFold
   cross_val.StratifiedKFold
   cross_val.LeaveOneLabelOut
   cross_val.LeavePLabelOut


Grid Search
===========

.. autosummary::
   :toctree: generated/
   :template: class.rst

   grid_search.GridSearchCV


.. _feature_selection_ref:


Feature Selection
=================

.. autosummary::
   :toctree: generated/
   :template: class.rst

   feature_selection.rfe.RFE
   feature_selection.rfe.RFECV


.. _feature_extraction_ref:


Feature Extraction
==================

.. autosummary::
   :toctree: generated/
   :template: function.rst

   feature_extraction.image.img_to_graph
   feature_extraction.image.grid_to_graph

.. autosummary::
   :toctree: generated/
   :template: class.rst

   feature_extraction.text.RomanPreprocessor
   feature_extraction.text.WordNGramAnalyzer
   feature_extraction.text.CharNGramAnalyzer
   feature_extraction.text.CountVectorizer
   feature_extraction.text.TfidfTransformer
   feature_extraction.text.Vectorizer

For sparse data
---------------

.. autosummary::
   :toctree: generated/
   :template: class.rst

   feature_extraction.text.sparse.TfidfTransformer
   feature_extraction.text.sparse.CountVectorizer
   feature_extraction.text.sparse.Vectorizer


Pipeline
========

.. autosummary::
   :toctree: generated/
   :template: class.rst

   pipeline.Pipeline


Partial Least Squares
=====================

.. autosummary::
   :toctree: generated/
   :template: class.rst

   pls.PLSRegression
   pls.PLSCanonical
   pls.CCA
   pls.PLSSVD

