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
-----------------

.. autosummary::
   :toctree: generated/
   :template: class.rst

   svm.sparse.SVC
   svm.sparse.NuSVC
   svm.sparse.SVR
   svm.sparse.NuSVR
   svm.sparse.OneClassSVM
   svm.sparse.LinearSVC

Logistic Regression
===================

.. autosummary::
   :toctree: generated/
   :template: class.rst

   logistic.LogisticRegression


Generalized Linear Models
=========================

.. autosummary::
   :toctree: generated/
   :template: class.rst

   glm.LinearRegression
   glm.Ridge
   glm.Lasso
   glm.LassoCV
   glm.ElasticNet
   glm.ElasticNetCV
   glm.LARS
   glm.LassoLARS


.. autosummary::
   :toctree: generated/
   :template: function.rst

   glm.lars_path


For sparse data
---------------

.. autosummary::
   :toctree: generated/
   :template: class.rst

   glm.sparse.Lasso
   glm.sparse.ElasticNet
        

Bayesian Regression
===================

.. autosummary::
   :toctree: generated/
   :template: class.rst

   glm.BayesianRidge
   glm.ARDRegression   

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

   neighbors.Neighbors
   ball_tree.BallTree

.. autosummary::
   :toctree: generated/
   :template: function.rst

   ball_tree.knn_brute

Gaussian Mixture Models
=======================

.. autosummary::
   :toctree: generated/
   :template: class.rst

   gmm.GMM


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


Covariance estimators
=====================


.. autosummary::
   :toctree: generated/
   :template: class.rst

   covariance.Covariance
   covariance.ShrunkCovariance
   covariance.LedoitWolf


Cross-validation
===================

.. autosummary::
   :toctree: generated/
   :template: class.rst
   
   cross_val.LeaveOneOut
   cross_val.LeavePOut
   cross_val.KFold
   cross_val.StratifiedKFold
   cross_val.LeaveOneLabelOut
   cross_val.LeavePLabelOut

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

.. autosummary::
   :toctree: generated/
   :template: class.rst

   feature_extraction.text.WordNGramAnalyzer
   feature_extraction.text.CharNGramAnalyzer
   feature_extraction.text.TermCountVectorizer
   feature_extraction.text.TfidfTransformer
   feature_extraction.text.TfidfVectorizer
   feature_extraction.text.SparseHashingVectorizer 

Pipeline
========

.. autosummary::
   :toctree: generated/
   :template: class.rst

   pipeline.Pipeline
