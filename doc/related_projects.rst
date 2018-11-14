.. _related_projects:

=====================================
Related Projects
=====================================

Projects implementing the scikit-learn estimator API are encouraged to use
the `scikit-learn-contrib template <https://github.com/scikit-learn-contrib/project-template>`_
which facilitates best practices for testing and documenting estimators.
The `scikit-learn-contrib GitHub organisation <https://github.com/scikit-learn-contrib/scikit-learn-contrib>`_
also accepts high-quality contributions of repositories conforming to this
template.

Below is a list of sister-projects, extensions and domain specific packages.

Interoperability and framework enhancements
-------------------------------------------

These tools adapt scikit-learn for use with other technologies or otherwise
enhance the functionality of scikit-learn's estimators.

**Data formats**

- `sklearn_pandas <https://github.com/paulgb/sklearn-pandas/>`_ bridge for
  scikit-learn pipelines and pandas data frame with dedicated transformers.
  
- `sklearn_xarray <https://github.com/phausamann/sklearn-xarray/>`_ provides
  compatibility of scikit-learn estimators with xarray data structures.

**Auto-ML**

- `auto_ml <https://github.com/ClimbsRocks/auto_ml/>`_
  Automated machine learning for production and analytics, built on scikit-learn
  and related projects. Trains a pipeline wth all the standard machine learning 
  steps. Tuned for prediction speed and ease of transfer to production environments. 

- `auto-sklearn <https://github.com/automl/auto-sklearn/>`_
  An automated machine learning toolkit and a drop-in replacement for a
  scikit-learn estimator

- `TPOT <https://github.com/rhiever/tpot>`_
  An automated machine learning toolkit that optimizes a series of scikit-learn
  operators to design a machine learning pipeline, including data and feature
  preprocessors as well as the estimators. Works as a drop-in replacement for a
  scikit-learn estimator.

- `scikit-optimize <https://scikit-optimize.github.io/>`_
  A library to minimize (very) expensive and noisy black-box functions. It
  implements several methods for sequential model-based optimization, and
  includes a replacement for ``GridSearchCV`` or ``RandomizedSearchCV`` to do
  cross-validated parameter search using any of these strategies.

**Experimentation frameworks**

- `REP <https://github.com/yandex/REP>`_ Environment for conducting data-driven
  research in a consistent and reproducible way

- `ML Frontend <https://github.com/jeff1evesque/machine-learning>`_ provides
  dataset management and SVM fitting/prediction through
  `web-based <https://github.com/jeff1evesque/machine-learning#web-interface>`_
  and `programmatic <https://github.com/jeff1evesque/machine-learning#programmatic-interface>`_
  interfaces.

- `Scikit-Learn Laboratory
  <https://skll.readthedocs.io/en/latest/index.html>`_  A command-line
  wrapper around scikit-learn that makes it easy to run machine learning
  experiments with multiple learners and large feature sets.

- `Xcessiv <https://github.com/reiinakano/xcessiv>`_ is a notebook-like
  application for quick, scalable, and automated hyperparameter tuning
  and stacked ensembling. Provides a framework for keeping track of 
  model-hyperparameter combinations.

**Model inspection and visualisation**

- `eli5 <https://github.com/TeamHG-Memex/eli5/>`_ A library for
  debugging/inspecting machine learning models and explaining their
  predictions.

- `mlxtend <https://github.com/rasbt/mlxtend>`_ Includes model visualization
  utilities.

- `scikit-plot <https://github.com/reiinakano/scikit-plot>`_ A visualization library
  for quick and easy generation of common plots in data analysis and machine learning.

- `yellowbrick <https://github.com/DistrictDataLabs/yellowbrick>`_ A suite of
  custom matplotlib visualizers for scikit-learn estimators to support visual feature
  analysis, model selection, evaluation, and diagnostics.


**Model export for production**

- `sklearn-pmml <https://github.com/alex-pirozhenko/sklearn-pmml>`_
  Serialization of (some) scikit-learn estimators into PMML.

- `sklearn2pmml <https://github.com/jpmml/sklearn2pmml>`_
  Serialization of a wide variety of scikit-learn estimators and transformers
  into PMML with the help of `JPMML-SkLearn <https://github.com/jpmml/jpmml-sklearn>`_
  library.

- `sklearn-porter <https://github.com/nok/sklearn-porter>`_
  Transpile trained scikit-learn models to C, Java, Javascript and others.

- `sklearn-compiledtrees <https://github.com/ajtulloch/sklearn-compiledtrees/>`_
  Generate a C++ implementation of the predict function for decision trees (and
  ensembles) trained by sklearn. Useful for latency-sensitive production
  environments.


Other estimators and tasks
--------------------------

Not everything belongs or is mature enough for the central scikit-learn
project. The following are projects providing interfaces similar to
scikit-learn for additional learning algorithms, infrastructures
and tasks.

**Structured learning**

- `Seqlearn <https://github.com/larsmans/seqlearn>`_  Sequence classification
  using HMMs or structured perceptron.

- `HMMLearn <https://github.com/hmmlearn/hmmlearn>`_ Implementation of hidden
  markov models that was previously part of scikit-learn.

- `PyStruct <https://pystruct.github.io>`_ General conditional random fields
  and structured prediction.

- `pomegranate <https://github.com/jmschrei/pomegranate>`_ Probabilistic modelling
  for Python, with an emphasis on hidden Markov models.

- `sklearn-crfsuite <https://github.com/TeamHG-Memex/sklearn-crfsuite>`_
  Linear-chain conditional random fields
  (`CRFsuite <http://www.chokkan.org/software/crfsuite/>`_ wrapper with
  sklearn-like API).

**Deep neural networks etc.**

- `pylearn2 <http://deeplearning.net/software/pylearn2/>`_ A deep learning and
  neural network library build on theano with scikit-learn like interface.

- `sklearn_theano <https://sklearn-theano.github.io/>`_ scikit-learn compatible
  estimators, transformers, and datasets which use Theano internally

- `nolearn <https://github.com/dnouri/nolearn>`_ A number of wrappers and
  abstractions around existing neural network libraries

- `keras <https://github.com/fchollet/keras>`_ Deep Learning library capable of
  running on top of either TensorFlow or Theano.

- `lasagne <https://github.com/Lasagne/Lasagne>`_ A lightweight library to
  build and train neural networks in Theano.
  
- `skorch <https://github.com/dnouri/skorch>`_ A scikit-learn compatible 
  neural network library that wraps PyTorch.

**Broad scope**

- `mlxtend <https://github.com/rasbt/mlxtend>`_ Includes a number of additional
  estimators as well as model visualization utilities.

- `sparkit-learn <https://github.com/lensacom/sparkit-learn>`_ Scikit-learn
  API and functionality for PySpark's distributed modelling.

**Other regression and classification**

- `xgboost <https://github.com/dmlc/xgboost>`_ Optimised gradient boosted decision
  tree library.

- `ML-Ensemble <https://mlens.readthedocs.io/>`_ Generalized
  ensemble learning (stacking, blending, subsemble, deep ensembles,
  etc.).

- `lightning <https://github.com/scikit-learn-contrib/lightning>`_ Fast
  state-of-the-art linear model solvers (SDCA, AdaGrad, SVRG, SAG, etc...).

- `py-earth <https://github.com/scikit-learn-contrib/py-earth>`_ Multivariate
  adaptive regression splines

- `Kernel Regression <https://github.com/jmetzen/kernel_regression>`_
  Implementation of Nadaraya-Watson kernel regression with automatic bandwidth
  selection

- `gplearn <https://github.com/trevorstephens/gplearn>`_ Genetic Programming
  for symbolic regression tasks.

- `multiisotonic <https://github.com/alexfields/multiisotonic>`_ Isotonic
  regression on multidimensional features.

- `scikit-multilearn <https://scikit.ml>`_ Multi-label classification with 
  focus on label space manipulation.

- `seglearn <https://github.com/dmbee/seglearn>`_ Time series and sequence 
  learning using sliding window segmentation.

**Decomposition and clustering**

- `lda <https://github.com/ariddell/lda/>`_: Fast implementation of latent
  Dirichlet allocation in Cython which uses `Gibbs sampling
  <https://en.wikipedia.org/wiki/Gibbs_sampling>`_ to sample from the true
  posterior distribution. (scikit-learn's
  :class:`sklearn.decomposition.LatentDirichletAllocation` implementation uses
  `variational inference
  <https://en.wikipedia.org/wiki/Variational_Bayesian_methods>`_ to sample from
  a tractable approximation of a topic model's posterior distribution.)

- `Sparse Filtering <https://github.com/jmetzen/sparse-filtering>`_
  Unsupervised feature learning based on sparse-filtering

- `kmodes <https://github.com/nicodv/kmodes>`_ k-modes clustering algorithm for
  categorical data, and several of its variations.

- `hdbscan <https://github.com/scikit-learn-contrib/hdbscan>`_ HDBSCAN and Robust Single
  Linkage clustering algorithms for robust variable density clustering.

- `spherecluster <https://github.com/clara-labs/spherecluster>`_ Spherical
  K-means and mixture of von Mises Fisher clustering routines for data on the
  unit hypersphere.

**Pre-processing**

- `categorical-encoding
  <https://github.com/scikit-learn-contrib/categorical-encoding>`_ A
  library of sklearn compatible categorical variable encoders.

- `imbalanced-learn
  <https://github.com/scikit-learn-contrib/imbalanced-learn>`_ Various
  methods to under- and over-sample datasets.

Statistical learning with Python
--------------------------------
Other packages useful for data analysis and machine learning.

- `Pandas <https://pandas.pydata.org/>`_ Tools for working with heterogeneous and
  columnar data, relational queries, time series and basic statistics.

- `theano <http://deeplearning.net/software/theano/>`_ A CPU/GPU array
  processing framework geared towards deep learning research.

- `statsmodels <https://www.statsmodels.org>`_ Estimating and analysing
  statistical models. More focused on statistical tests and less on prediction
  than scikit-learn.

- `PyMC <https://pymc-devs.github.io/pymc/>`_ Bayesian statistical models and
  fitting algorithms.

- `Sacred <https://github.com/IDSIA/Sacred>`_ Tool to help you configure,
  organize, log and reproduce experiments

- `Seaborn <https://stanford.edu/~mwaskom/software/seaborn/>`_ Visualization library based on
  matplotlib. It provides a high-level interface for drawing attractive statistical graphics.

- `Deep Learning <http://deeplearning.net/software_links/>`_ A curated list of deep learning
  software libraries.


Domain specific packages
~~~~~~~~~~~~~~~~~~~~~~~~

- `scikit-image <https://scikit-image.org/>`_ Image processing and computer
  vision in python.

- `Natural language toolkit (nltk) <https://www.nltk.org/>`_ Natural language
  processing and some machine learning.

- `gensim <https://radimrehurek.com/gensim/>`_  A library for topic modelling,
  document indexing and similarity retrieval

- `NiLearn <https://nilearn.github.io/>`_ Machine learning for neuro-imaging.

- `AstroML <https://www.astroml.org/>`_  Machine learning for astronomy.

- `MSMBuilder <http://msmbuilder.org/>`_  Machine learning for protein
  conformational dynamics time series.

- `scikit-surprise <https://surpriselib.com/>`_ A scikit for building and
  evaluating recommender systems.

Snippets and tidbits
---------------------

The `wiki <https://github.com/scikit-learn/scikit-learn/wiki/Third-party-projects-and-code-snippets>`_ has more!
