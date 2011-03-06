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
expect a numpy array as input ``X``.  The expected shape of X is
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


 :Features in the extended Iris dataset:

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


How to build "good" feature extraction strategy
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The good rule of thumb is to imagine a human being performing the
task the machine is trying to accomplish using only the numerical
features provided to the machine.

Usually it helps if and only if two samples **judged similar in
real life** by the human being are **close according to some
similarity metric in feature space**.

In other words, the feature extraction strategy must somehow preserve
the intuitive topology of the sample set.


Supervised Learning: ``model.fit(X, y)``
----------------------------------------

 - Principles

 - Sample algorithms in ``scikit-learn``

 - Real life applications


Unsupervised Learning: ``model.fit(X)``
---------------------------------------

 - Principles

 - Sample algorithms in ``scikit-learn``

 - Real life applications



Training set, test sets and overfitting
---------------------------------------

