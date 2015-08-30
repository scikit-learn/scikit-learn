.. _faq:

===========================
Frequently Asked Questions
===========================

Here we try to give some answers to questions that regularly pop up on the mailing list.

What is the project name (a lot of people get it wrong)?
--------------------------------------------------------
scikit-learn, but not scikit or SciKit nor sci-kit learn. Also not scikits.learn or scikits-learn, which where previously used.

How do you pronounce the project name?
------------------------------------------
sy-kit learn. sci stands for science!

Why scikit?
------------
There are multiple scikits, which are scientific toolboxes build around SciPy.
You can find a list at `<https://scikits.appspot.com/scikits>`_.
Apart from scikit-learn, another popular one is `scikit-image <http://scikit-image.org/>`_.

How can I contribute to scikit-learn?
-----------------------------------------
See :ref:`contributing`. Before wanting to add a new algorithm, which is
usually a major and lengthy undertaking, it is recommended to start with :ref:`known
issues <easy_issues>`.


How can I create a bunch object?
------------------------------------------------

Don't make a bunch object! They are not part of the scikit-learn API. Bunch
objects are just a way to package some numpy arrays. As a scikit-learn user you
only ever need numpy arrays to feed your model with data.

For instance to train a classifier, all you need is a 2D array ``X`` for the
input variables and a 1D array ``y`` for the target variables. The array ``X``
holds the features as columns and samples as rows . The array ``y`` contains
integer values to encode the class membership of each sample in ``X``.

To load data as numpy arrays you can use different libraries depending on the
original data format:

* `numpy.loadtxt
  <http://docs.scipy.org/doc/numpy/reference/generated/numpy.loadtxt.html>`_ to
  load text files (such as CSV) assuming that all the columns have an
  homogeneous data type (e.g. all numeric values).

* `scipy.io <http://docs.scipy.org/doc/scipy/reference/io.html>`_ for common
  binary formats often used in scientific computing context.

* `scipy.misc.imread <http://docs.scipy.org/doc/scipy/reference/generated/scipy.
  misc.imread.html#scipy.misc.imread>`_ (requires the `Pillow
  <https://pypi.python.org/pypi/Pillow>`_ package) to load pixel intensities
  data from various image file formats.

* `pandas.io <http://pandas.pydata.org/pandas-docs/stable/io.html>`_ to load
  heterogeneously typed data from various file formats and database protocols
  that can slice and dice before conversion to numerical features in a numpy
  array.

Note: if you manage your own numerical data it is recommended to use an
optimized file format such as HDF5 to reduce data load times. Various libraries
such as H5Py, PyTables and pandas provides a Python interface for reading and
writing data in that format.

Can I add this new algorithm that I (or someone else) just published?
-------------------------------------------------------------------------
No. As a rule we only add well-established algorithms. A rule of thumb is at least
3 years since publications, 200+ citations and wide use and usefullness. A
technique that provides a clear-cut improvement (e.g. an enhanced  data
structure or efficient approximation) on a widely-used method will also be
considered for inclusion.
Your implementation doesn't need to be in scikit-learn to be used together
with scikit-learn tools, though. Implement your favorite algorithm
in a scikit-learn compatible way, upload it to github and we will list
it under :ref:`related_projects`.
Also see selectiveness_.


Can I add this classical algorithm from the 80s?
---------------------------------------------------
Depends. If there is a common usecase within the scope of scikit-learn, such
as classification, regression or clustering, where it outperforms methods
that are already implemented in scikit-learn, we will consider it.

.. _selectiveness:

Why are you so selective on what algorithms you include in scikit-learn?
------------------------------------------------------------------------
Code is maintenance cost, and we need to balance the amount of
code we have with the size of the team (and add to this the fact that
complexity scales non linearly with the number of features).
The package relies on core developers using their free time to
fix bugs, maintain code and review contributions.
Any algorithm that is added needs future attention by the developers,
at which point the original author might long have lost interest.
Also see `this thread on the mailing list
<http://sourceforge.net/p/scikit-learn/mailman/scikit-learn-general/thread/CAAkaFLWcBG%2BgtsFQzpTLfZoCsHMDv9UG5WaqT0LwUApte0TVzg%40mail.gmail.com/#msg33104380>`_.

Why did you remove HMMs from scikit-learn?
--------------------------------------------
See :ref:`adding_graphical_models`.

.. _adding_graphical_models:

Will you add graphical models or sequence prediction to scikit-learn?
------------------------------------------------------------------------
Not in the foreseeable future.
scikit-learn tries to provide a unified API for the basic tasks in machine
learning, with pipelines and meta-algorithms like grid search to tie
everything together. The required concepts, APIs, algorithms and
expertise required for structured learning are different from what
scikit-learn has to offer. If we started doing arbitrary structured
learning, we'd need to redesign the whole package and the project
would likely collapse under its own weight.

There are two project with API similar to scikit-learn that
do structured prediction:

* `pystruct <http://pystruct.github.io/>`_ handles general structured
  learning (focuses on SSVMs on arbitrary graph structures with
  approximate inference; defines the notion of sample as an instance of
  the graph structure)

* `seqlearn <http://larsmans.github.io/seqlearn/>`_ handles sequences only (focuses on
  exact inference; has HMMs, but mostly for the sake of completeness;
  treats a feature vector as a sample and uses an offset encoding for
  the dependencies between feature vectors)

Will you add GPU support?
--------------------------
No, or at least not in the near future. The main reason is that GPU support will introduce many software dependencies and introduce platform specific issues.
scikit-learn is designed to be easy to install on a wide variety of platforms.
Outside of neural networks, GPUs don't play a large role in machine learning today, and much larger gains in speed can often be achieved by a careful choice of algorithms.

Do you support PyPy?
--------------------

In case you didn't know, `PyPy <http://pypy.org/>`_ is the new, fast,
just-in-time compiling Python implementation. We don't support it.
When the `NumPy support <http://buildbot.pypy.org/numpy-status/latest.html>`_
in PyPy is complete or near-complete, and SciPy is ported over as well,
we can start thinking of a port.
We use too much of NumPy to work with a partial implementation.

How do I deal with string data (or trees, graphs...)?
-----------------------------------------------------

scikit-learn estimators assume you'll feed them real-valued feature vectors.
This assumption is hard-coded in pretty much all of the library.
However, you can feed non-numerical inputs to estimators in several ways.

If you have text documents, you can use a term frequency features; see
:ref:`text_feature_extraction` for the built-in *text vectorizers*.
For more general feature extraction from any kind of data, see
:ref:`dict_feature_extraction` and :ref:`feature_hashing`.

Another common case is when you have non-numerical data and a custom distance
(or similarity) metric on these data. Examples include strings with edit
distance (aka. Levenshtein distance; e.g., DNA or RNA sequences). These can be
encoded as numbers, but doing so is painful and error-prone. Working with
distance metrics on arbitrary data can be done in two ways.

Firstly, many estimators take precomputed distance/similarity matrices, so if
the dataset is not too large, you can compute distances for all pairs of inputs.
If the dataset is large, you can use feature vectors with only one "feature",
which is an index into a separate data structure, and supply a custom metric
function that looks up the actual data in this data structure. E.g., to use
DBSCAN with Levenshtein distances::

    >>> from leven import levenshtein       # doctest: +SKIP
    >>> import numpy as np
    >>> from sklearn.cluster import dbscan
    >>> data = ["ACCTCCTAGAAG", "ACCTACTAGAAGTT", "GAATATTAGGCCGA"]
    >>> def lev_metric(x, y):
    ...     i, j = int(x[0]), int(y[0])     # extract indices
    ...     return levenshtein(data[i], data[j])
    ...
    >>> X = np.arange(len(data)).reshape(-1, 1)
    >>> X
    array([[0],
           [1],
           [2]])
    >>> dbscan(X, metric=lev_metric, eps=5, min_samples=2)  # doctest: +SKIP
    ([0, 1], array([ 0,  0, -1]))

(This uses the third-party edit distance package ``leven``.)

Similar tricks can be used, with some care, for tree kernels, graph kernels,
etc.
