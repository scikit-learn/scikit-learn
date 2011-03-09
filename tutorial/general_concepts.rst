Machine Learning 101
====================

Machine Learning is about building **programs with tunable parameters**
(typically an array of floating point values) that are adjusted
automatically so as to improve its behavior by **adapting to
previously seen data**.

Machine Learning can be considered a **subfield of Artificial
Intelligence** since those algorithms can be seen as building blocks
to make computer learn to behave more intelligently by somehow
**generalizing** rather that just storing and retrieving data items
like a database system would do.

The following will introduce the main concepts used to qualify
machine learning algorithms, show how those concepts match with the
``scikit-learn`` API and give example applications.


Features and feature extraction
-------------------------------

Most machine learning algorithms implemented in ``scikit-learn``
expect a numpy array as input ``X``.  The expected shape of ``X`` is
``(n_samples, n_features)``.

:``n_samples``:

  The number of samples: each sample is an item to process (e.g.
  classifiy). A sample can be a document, a picture, a sound, a
  video, a row in database or CSV file, or whatever you can
  describe with a fixed set of quantitative traits.

:``n_features``:

  The number of features or distinct traits that can be used to
  describe each item in a quantitative manner.


The number of features must be fixed in advance. However it can be
very high dimensional (e.g. millions of features) with most of them
being zeros for a given sample. In this case we use ``scipy.sparse``
matrices instead of ``numpy`` arrays so has to make the data fit
in memory.


A simple example: the iris dataset
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The machine learning community often uses a simple flowers database
were each row in the database (or CSV file) is a set of measurements
of an individual iris flower.

.. figure:: images/Virginia_Iris.png
   :scale: 100 %
   :align: center
   :alt: Photo of Iris Virginia

   Iris Virginia (source: Wikipedia)


Each sample in this dataset is described by 4 features and can
belong to one of the target classes:

 :Features in the Iris dataset:

   0. sepal length in cm
   1. sepal width in cm
   2. petal length in cm
   3. petal width in cm

 :Target classes to predict:

   0. Iris Setosa
   1. Iris Versicolour
   2. Iris Virginica


``scikit-learn`` embeds a copy of the iris CSV file along with a
helper function to load it into numpy arrays::

  >>> from scikits.learn.datasets import load_iris
  >>> iris = load_iris()

.. note::

  To be able to copy and paste examples without taking care of the leading
  ``>>>`` and ``...`` prompt signs, enable the ipython doctest mode with:
  ``%doctest_mode``

The features of each sample flower is stored in the ``data`` attribute
of the dataset::

  >>> n_samples, n_features = iris.data.shape

  >>> n_samples
  150

  >>> n_features
  4

  >>> iris.data[0]
  array([ 5.1,  3.5,  1.4,  0.2])


The information about the class of each sample is stored in the
``target`` attribute of the dataset::

  >>> len(iris.target) == n_samples
  True

  >>> iris.target
  array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
         0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
         1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
         1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
         2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
         2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2])

The names of the classes is stored in the last attribute, namely
``target_names``::

  >>> list(iris.target_names)
  ['setosa', 'versicolor', 'virginica']


Handling categorical features
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Sometimes people describe samples with categorical descriptors that
have no obvious numerical representation. For instance assume that
each flower is further described by a color name among a fixed list
of color names::

  color in ['purple', 'blue', 'red']

The simple way to turn this categorical feature into numerical
features suitable for machine learning is to create new features
for each distinct color name that can be valued to ``1.0`` if the
category is matching or ``0.0`` if not.

The enriched iris feature set would hence be in this case:

  0. sepal length in cm
  1. sepal width in cm
  2. petal length in cm
  3. petal width in cm
  4. color#purple (1.0 or 0.0)
  5. color#blue (1.0 or 0.0)
  6. color#red (1.0 or 0.0)


Extracting features from unstructured data
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The previous example deals with features that are readily available
in a structured datasets with rows and columns of numerical or
categorical values.

However, **most of the produced data is not readily available in a
structured representation** such as SQL, CSV, XML, JSON or RDF.

Here is an overview of strategies to turn unstructed data items
into arrays of numerical features.


  :Text documents:

    Count the frequency of each word or pair of consecutive words
    in each document. This approach is called the **Bag of Words**.

    Note: we include other files formats such as HTML and PDF in
    this category: an ad-hoc preprocessing step is required to
    extract the plain text in UTF-8 encoding for instance.


  :Images:

    - Rescale the picture to a fixed size and **take all the raw
      pixels values** (with or without luminosity normalization)

    - Take some transformation of the signal (gradients in each
      pixel, wavelets transforms...)

    - Compute the Euclidean, Manhattan or cosine **similarities of
      the sample to a set reference prototype images** aranged in a
      code book.  The code book may have been previously extracted
      on the same dataset using an unsupervised learning algorithms
      on the raw pixel signal.

      Each feature value is the distance to one element of the code
      book.

    - Perform **local feature extraction**: split the picture into
      small regions and perform feature extraction locally in each
      area.

      Then combine all the feature of the individual areas into a
      single array.

  :Sounds:

    Same strategy as for images with in a 1D space instead of 2D


Practical implementations of such feature extraction strategies
will be presented in the last sections of this tutorial.


How to devise a "good" feature extraction strategy
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The feature extraction strategy both depends on the task we are
trying to perform and the nature of the collected data. Therefore
there is no formal rule to define which strategy is the best.

A good rule of thumb is to imagine a human-being performing the
task the machine is trying to accomplish using only the numerical
features provided to the machine.

Usually the feature extraction is useful if and only if two samples
**judged similar in real life** by the human-being are **close
according to some similarity metric of the feature space**.

In other words, the feature extraction strategy must somehow preserve
the intuitive topology of the sample set.


Supervised Learning: ``model.fit(X, y)``
----------------------------------------

.. figure:: images/supervised.png
   :scale: 75 %
   :align: center
   :alt: Flow diagram for supervised learning

   Supervised Learning overview

A supervised learning algorithm makes the distinction between the
raw observed data ``X`` with shape ``(n_samples, n_features)`` and
some label given to the model while training by some teacher. In
``scikit-learn`` this array is often noted ``y`` and has generally
the shape ``(n_samples,)``.

After training, the fitted model does no longer expect the ``y``
as an input: it will try to predict the most likely labels ``y_new``
for new a set of samples ``X_new``.

Depending on the nature of the target ``y``, supervised learning
can be given different names:

  - If ``y`` has values in a fixed set of categorical outcomes
    (represented by integers) the task to predict ``y`` is called
    classification.

  - If ``y`` has floating point values (e.g. to represent a price,
    a temperature, a size...), the task to predict ``y`` is called
    regression.


Classification
~~~~~~~~~~~~~~


A first classifier example with ``scikit-learn``
++++++++++++++++++++++++++++++++++++++++++++++++

In the iris dataset example, suppose we are assigned the task to
guess the class of an individual flower given the measurements of
petals and sepals. This is a classification task, hence we have::

  >>> X, y = iris.data, iris.target

Once the data has this format it is trivial to train a classifier,
for instance a support vector machine with a linear kernel (or lack
of thereof)::

  >>> from scikits.learn.svm import LinearSVC
  >>> clf = LinearSVC()

``clf`` is a statistical model that has parameters that control the
learning algorithm (those parameters are sometimes called the
hyper-parameters). Those hyperparameters can be supplied by the
user in the constructore of the model. We will explain later choose
a good combination either using simple empirical rules or data
driven selection::

  >>> clf
  LinearSVC(loss='l2', C=1.0, intercept_scaling=1, fit_intercept=True,
       eps=0.0001, penalty='l2', multi_class=False, dual=True)

By default the real model parameters are not initialized. They will be
automatically be tuned from the data by calling the ``fit`` method::

  >>> clf = clf.fit(X, y)

  >>> clf.coef_
  array([[ 0.18423474,  0.45122764, -0.80794654, -0.45071379],
         [ 0.04864394, -0.88914385,  0.40540293, -0.93720122],
         [-0.85086062, -0.98671553,  1.38098573,  1.8653574 ]])

  >>> clf.intercept_
  array([ 0.10956015,  1.6738296 , -1.70973044])

Once the model is trained, it can be used to predict the most likely outcome on
unseen data. For instance let us define a list of a simple sample that looks
like the first sample of the iris dataset::

  >>> X_new = [[ 5.0,  3.6,  1.3,  0.25]]

  >>> clf.predict(X_new)
  array([0], dtype=int32)

The outcome is ``0`` which the id of the first iris class namely
'setosa'.

The following figure places the location of the fit and predict
calls on the previous flow diagram. The ``vec`` object is a vectorizer
used for feature extractor that is not used in the case of the iris
data which already comes as vectors of features:

.. figure:: images/supervised_scikit_learn.png
   :scale: 75 %
   :align: center
   :alt: Flow diagram for supervised learning with scikit-learn

   Supervised Learning with scikit-learn



Some ``scikit-learn`` classifiers can further predicts probabilities
of the outcome.  This is the case of logistic regression models::

  >>> from scikits.learn.linear_model import LogisticRegression
  >>> clf2 = LogisticRegression().fit(X, y)
  >>> clf2
  LogisticRegression(C=1.0, intercept_scaling=1, fit_intercept=True, eps=0.0001,
            penalty='l2', dual=False)

  >>> clf2.predict_proba(X_new)
  array([[  9.07512928e-01,   9.24770379e-02,   1.00343962e-05]])

This means that the model estimates that the sample in ``X_new`` has:

  - 90% likelyhood to be belong to the 'setosa' class

  - 9% likelyhood to be belong to the 'versicolor' class

  - 1% likelyhood to be belong to the 'virginica' class

Of course the ``predict`` method that output the label id of the
most likely outcome is also available::

  >>> clf2.predict(X_new)
  array([0], dtype=int32)


Notable implementations of classifiers
++++++++++++++++++++++++++++++++++++++

:``scikits.learn.linear_model.LogisticRegression``:

  Regularized Logistic Regression based on ``liblinear``

:``scikits.learn.svm.LinearSVC``:

  Support Vector Machines without kernels based on ``liblinear``

:``scikits.learn.svm.SVC``:

  Support Vector Machines with kernels based on ``libsvm``

:``scikits.learn.linear_model.SGDClassifier``:

  Regularized linear models (SVM or logistic regression) using a Stochastic
  Gradient Descent algorithm written in ``Cython``

:``scikits.learn.neighbors.NeighborsClassifier``:

  k-Nearest Neighbors classifier based on the ball tree datastructure for low
  dimensional data and brute force search for high dimensional data


Sample application of classifiers
+++++++++++++++++++++++++++++++++

The following table gives example of application of such algorithms for
some common engineering tasks:

============================================ =================================
Task                                         Predicted outcomes
============================================ =================================
E-mail classification                        Spam, normal, priority mail
-------------------------------------------- ---------------------------------
Language identification in text documents    en, es, de, fr, ja, zh, ar, ru...
-------------------------------------------- ---------------------------------
News articles categorization                 Business, technology, sports...
-------------------------------------------- ---------------------------------
Sentiment Analysis in customer feedback      Negative, neutral, positive
-------------------------------------------- ---------------------------------
Face verification in pictures                Same / different persons
-------------------------------------------- ---------------------------------
Speaker verification on voice recordings     Same / different persons
============================================ =================================


Regression
~~~~~~~~~~

TODO


Unsupervised Learning: ``model.fit(X)``
---------------------------------------

.. figure:: images/unsupervised.png
   :scale: 75 %
   :align: center
   :alt: Flow diagram for unsupervised learning

   Unsupervised Learning overview

An unsupervised learning algorithm only uses a single set of
observations ``X`` with shape ``(n_samples, n_features)`` and does
not use any kind of labels.

An unsupervised learning model will try to fit its parameters so
as to best summarize regularities found in the data.

The following introduces the main variants of unsupervised learning
algorithms.


Dimensionality Reduction and visualization
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Dimensionality reduction the task to **derive a set of new artificial
features that is smaller than the original feature set while retaining
most of the variance of the original data**.


Normalization and visualization with PCA
++++++++++++++++++++++++++++++++++++++++

The most common technique for dimensionality reduction is called
**Principal Component Analysis**.

PCA can be done using linear combinations of the original features
using a truncated Singular Value Decomposition of the matrix ``X``
so as to project the data onto a base of the top singular vectors.

If the number of retained components is 2 or 3, PCA can be used to
visualize the dataset::


  >>> from scikits.learn.pca import PCA
  >>> pca = PCA(n_components=2, whiten=True).fit(X)

Once fitted, the ``pca`` model exposes the singular vectors as in the
``components_`` attribute::

  >>> pca.components_.T
  array([[ 0.17650757, -0.04015901,  0.41812992,  0.17516725],
         [-1.33840478, -1.48757227,  0.35831476,  0.15229463]])

  >>> pca.explained_variance_ratio_
  array([ 0.92461621,  0.05301557])

  >>> pca.explained_variance_ratio_.sum()
  0.97763177502480336

Let us project the iris dataset along those first 3 dimensions::

  >>> X_pca = pca.transform(X)

The dataset has been "normalized", which means that the data is now centered on
both components with unit variance::

  >>> X_pca.mean(axis=0)
  array([ -1.42478621e-15,   1.71936539e-15])

  >>> X_pca.std(axis=0)
  array([ 1.,  1.])

Furthermore the samples components do no longer carry any linear
correlation::

  >>> import numpy as np
  >>> np.corrcoef(X_pca.T)
  array([[  1.00000000e+00,   4.60742555e-16],
         [  4.60742555e-16,   1.00000000e+00]])


And visualize the dataset using ``pylab``, for instance by defining the
following utility function::

  >>> import pylab as pl
  >>> from itertools import cycle
  >>> def plot_2D(data, target, target_names):
  ...     colors = cycle('rgbcmykw')
  ...     target_ids = range(len(target_names))
  ...     pl.figure()
  ...     for i, c, label in zip(target_ids, colors, target_names):
  ...         pl.scatter(data[target == i, 0], data[target == i, 1],
  ...                    c=c, label=label)
  ...     pl.legend()
  ...

Calling ``plot_2D(X_pca, iris.target, iris.target_names)`` will
display the following:


.. figure:: images/iris_pca_2d.png
   :scale: 65 %
   :align: center
   :alt: 2D PCA projection of the iris dataset

   2D PCA projection of the iris dataset


.. note::

  The default implementation of PCA computes the SVD of the full
  data matrix which is not scalable when both ``n_samples`` and
  ``n_features`` are big (more that a few thousands).

  If you are interested in a number of components that is much
  smaller than both ``n_samples`` and ``n_features`` consider using
  ``scikits.learn.pca.RandomizedPCA`` instead.


Other applications of dimensionality reduction
++++++++++++++++++++++++++++++++++++++++++++++

Dimensionality Reduction is not just useful for visualization of
high dimensional datasets. It can also be used as a preprocessing
step (often called data normalization) to help speed up supervised
machine learning that are not computationally efficient with high
``n_features`` such as SVM classifiers with gaussian kernels for
instance or that do not work well with linearly correlated features.

.. note::

  ``scikit-learn`` also features an implementation of Independant
  Component Analysis (ICA) and work is under way to implement common
  manifold extraction strategies.


Clustering
~~~~~~~~~~

Clustering is the task of gathering samples into groups of similar
samples according to some predifined similarity or dissimilarity
measure (such as the Euclidean distance).

For instance let us reuse the output of the 2D PCA of the iris
dataset and try to find 3 groups of samples using the slimplest
clustering algorithm (KMeans)::

  >>> from scikits.learn.cluster import KMeans
  >>> from numpy.random import RandomState
  >>> rng = RandomState(42)

  >>> kmeans = KMeans(3, rng=rng).fit(X_pca)

  >>> kmeans.cluster_centers_
  array([[ 1.01505989, -0.70632886],
         [ 0.33475124,  0.89126382],
         [-1.287003  , -0.43512572]])


  >>> kmeans.labels_[:10]
  array([2, 2, 2, 2, 2, 2, 2, 2, 2, 2])

  >>> kmeans.labels_[-10:]
  array([0, 0, 1, 0, 0, 0, 1, 0, 0, 1])

We can plot the assigned cluster labels instead of the target names
with::

   plot_2D(X_pca, kmeans.labels_, ["c0", "c1", "c2"])


.. figure:: images/iris_pca_2d_kmeans.png
   :scale: 65 %
   :align: center
   :alt: KMeans cluster assignements on 2D PCA iris data

   KMeans cluster assignements on 2D PCA iris data


Notable implementations of clustering models
++++++++++++++++++++++++++++++++++++++++++++

The following are two well known clustering algorithms. Like most
unsupervised learning models in the scikit, they expect the data
to be cluster to have shape ``(n_samples, n_features)``:

:``scikits.learn.cluster.KMeans``:

  The simplest yet effective clustering algorithm. Need to be
  provided the number of clusters in advance and assume that the
  data is normalized as input (but use a PCA model as preprocessor).

:``scikits.learn.cluster.MeanShift``:

  Can find better looking clusters than KMeans but is not scalable
  to high number of samples.

Other clustering algorithms do not work with a data matrix with
shape ``(n_samples, n_features)`` but directly with a precomputed
affinity matrix with shape ``(n_samples, n_samples)``:

:``scikits.learn.cluster.AffinityPropagation``:

  Clustering algorithm based on message passing between data points.

:``scikits.learn.cluster.SpectralClustering``:

  KMeans applied to a projection of the normalized graph Laplacian:
  finds normalized graph cuts if the affinity matrix is interpreted
  as an adjacency matrix of a graph.


Hierarchical clustering is being implemented in a branch that is
likely to be merged into master before the release of ``scikit-learn``
0.8.


Applications of clustering
++++++++++++++++++++++++++

Here are some common applications of clustering algorithms:

 - Building customer profiles for market analysis
 - Grouping related web news (e.g. Google News) and websearch results
 - Grouping related stock quotes for investment portfolio management
 - Can be used as a preprocessing step for recommender systems
 - Can be used to build a code book of prototype samples for unsupervised
   feature extraction for supervised learning algorithms


Density estimation and outliers detection
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


TODO


Linearly separable data
-----------------------

- TODO: explain the concept using the 2D projected iris data

- Play with the interactive example from the ``examples`` folder of the
  ``scikit-learn`` distribution:

    % python $SKL_HOME/examples/applications/svm_gui.py

- Exercise: find a model that is able to solve the XOR problem using the GUI


Training set, test sets and overfitting
---------------------------------------

TODO


