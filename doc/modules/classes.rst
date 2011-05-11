===============
Class reference
===============


.. TODO: put some order here. Alphabetical ?


Support Vector Machines
=======================


.. automodule:: scikits.learn.svm
   :no-members:
   :no-inherited-members:

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

.. autosummary::

   :toctree: generated/
   :template: function.rst

   svm.l1_min_c


For sparse data
---------------

.. automodule:: scikits.learn.svm.sparse
   :no-members:
   :no-inherited-members:

.. currentmodule:: scikits.learn

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

.. automodule:: scikits.learn.linear_model
   :no-members:
   :no-inherited-members:

.. currentmodule:: scikits.learn

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
   linear_model.BayesianRidge
   linear_model.ARDRegression


.. autosummary::

   :toctree: generated/
   :template: function.rst

   linear_model.lasso_path
   linear_model.lars_path


For sparse data
---------------

.. automodule:: scikits.learn.linear_model.sparse
   :no-members:
   :no-inherited-members:

.. currentmodule:: scikits.learn

.. autosummary::

   :toctree: generated/
   :template: class.rst

   linear_model.sparse.Lasso
   linear_model.sparse.ElasticNet
   linear_model.sparse.SGDClassifier
   linear_model.sparse.SGDRegressor


Naive Bayes
===========

.. automodule:: scikits.learn.naive_bayes
   :no-members:
   :no-inherited-members:

.. currentmodule:: scikits.learn

.. autosummary::

   :toctree: generated/
   :template: class.rst

   naive_bayes.GNB


Nearest Neighbors
=================

.. automodule:: scikits.learn.neighbors
   :no-members:
   :no-inherited-members:

.. currentmodule:: scikits.learn

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

.. automodule:: scikits.learn.mixture
   :no-members:
   :no-inherited-members:

.. currentmodule:: scikits.learn

.. autosummary::

   :toctree: generated/
   :template: class.rst

   mixture.GMM


Hidden Markov Models
====================

.. automodule:: scikits.learn.hmm
   :no-members:
   :no-inherited-members:

.. currentmodule:: scikits.learn

.. autosummary::

   :toctree: generated/
   :template: class.rst

   hmm.GaussianHMM
   hmm.MultinomialHMM
   hmm.GMMHMM


Clustering
==========

.. automodule:: scikits.learn.cluster
   :no-members:
   :no-inherited-members:

.. currentmodule:: scikits.learn

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

.. automodule:: scikits.learn.metrics
   :no-members:
   :no-inherited-members:

.. currentmodule:: scikits.learn

.. autosummary::

   :toctree: generated/
   :template: function.rst

   metrics.euclidean_distances
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

Pairwise metrics
------------------

.. automodule:: scikits.learn.metrics.pairwise
   :no-members:
   :no-inherited-members:

.. currentmodule:: scikits.learn

.. autosummary::

   :toctree: generated/
   :template: function.rst

   metrics.pairwise.euclidean_distances
   metrics.pairwise.linear_kernel
   metrics.pairwise.polynomial_kernel
   metrics.pairwise.rbf_kernel


Covariance Estimators
=====================

.. automodule:: scikits.learn.covariance
   :no-members:
   :no-inherited-members:

.. currentmodule:: scikits.learn

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
   covariance.shrunk_covariance
   covariance.oas


Signal Decomposition
====================

.. automodule:: scikits.learn.decomposition
   :no-members:
   :no-inherited-members:

.. currentmodule:: scikits.learn

.. autosummary::

   :toctree: generated/
   :template: class.rst

   decomposition.PCA
   decomposition.ProbabilisticPCA
   decomposition.RandomizedPCA
   decomposition.KernelPCA
   decomposition.FastICA
   decomposition.NMF

.. autosummary::

   :toctree: generated/
   :template: function.rst

   decomposition.fastica


Linear Discriminant Analysis
============================

.. autosummary::

   :toctree: generated
   :template: class.rst

   lda.LDA


Cross Validation
================

.. automodule:: scikits.learn.cross_val
   :no-members:
   :no-inherited-members:

.. currentmodule:: scikits.learn

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

.. automodule:: scikits.learn.grid_search
   :no-members:
   :no-inherited-members:

.. currentmodule:: scikits.learn

.. autosummary::

   :toctree: generated/
   :template: class.rst

   grid_search.GridSearchCV


.. _feature_selection_ref:


Feature Selection
=================

.. automodule:: scikits.learn.feature_selection
   :no-members:
   :no-inherited-members:

.. currentmodule:: scikits.learn

.. autosummary::

   :toctree: generated/
   :template: class.rst

   feature_selection.rfe.RFE
   feature_selection.rfe.RFECV


.. _feature_extraction_ref:


Feature Extraction
==================

.. automodule:: scikits.learn.feature_extraction
   :no-members:
   :no-inherited-members:

.. currentmodule:: scikits.learn

From images
-------------

.. automodule:: scikits.learn.feature_extraction.image
   :no-members:
   :no-inherited-members:

.. currentmodule:: scikits.learn

.. autosummary::

   :toctree: generated/
   :template: function.rst

   feature_extraction.image.img_to_graph
   feature_extraction.image.grid_to_graph


From text
-------------

.. automodule:: scikits.learn.feature_extraction.text
   :no-members:
   :no-inherited-members:

.. currentmodule:: scikits.learn

.. autosummary::

   :toctree: generated/
   :template: class.rst

   feature_extraction.text.RomanPreprocessor
   feature_extraction.text.WordNGramAnalyzer
   feature_extraction.text.CharNGramAnalyzer
   feature_extraction.text.CountVectorizer
   feature_extraction.text.TfidfTransformer
   feature_extraction.text.Vectorizer


Pipeline
========

.. automodule:: scikits.learn.pipeline
   :no-members:
   :no-inherited-members:

.. currentmodule:: scikits.learn

.. autosummary::

   :toctree: generated/
   :template: class.rst

   pipeline.Pipeline


Partial Least Squares
=====================

.. automodule:: scikits.learn.pls
   :no-members:
   :no-inherited-members:

.. currentmodule:: scikits.learn

.. autosummary::

   :toctree: generated/
   :template: class.rst

   pls.PLSRegression
   pls.PLSCanonical
   pls.CCA
   pls.PLSSVD

