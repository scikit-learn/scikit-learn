.. _contributing:

============
Contributing
============

This project is a community effort, and everyone is welcome to
contribute.

The project is hosted on http://github.com/scikit-learn/scikit-learn


Submitting a bug report
=======================

In case you experience issues using the package, do not hesitate
to submit a ticket to the
`Bug Tracker <http://github.com/scikit-learn/scikit-learn/issues>`_.

You are also welcome to post there feature requests or links to pull-requests.


.. _git_repo:

Retrieving the latest code
==========================

You can check the latest sources with the command::

    git clone git://github.com/scikit-learn/scikit-learn.git

or if you have write privileges::

    git clone git@github.com:scikit-learn/scikit-learn.git

If you run the development version, it is cumbersome to re-install the
package each time you update the sources. It is thus preferred that
you add the scikit-directory to your ``PYTHONPATH`` and build the
extension in place::

    python setup.py build_ext --inplace

On Unix you can simply type ``make`` in the top-level folder to build
in-place and launch all the tests. Have a look at the ``Makefile`` for
additional utilities.


Contributing code
=================

.. note::

  To avoid duplicated work it is highly advised to contact the developers
  mailing list before starting work on a non-trivial feature.

  https://lists.sourceforge.net/lists/listinfo/scikit-learn-general

How to contribute
-----------------

The prefered way to contribute to Scikit-Learn is to fork the main
repository on
`github <http://github.com/scikit-learn/scikit-learn/>`__:

 1. `Create an account <https://github.com/signup/free>`_ on
    github if you don't have one already.

 2. Fork the `project repository
    <http://github.com/scikit-learn/scikit-learn>`__: click on the 'Fork'
    button, at the top, center of the page. This creates a copy of
    the code on the GitHub server where you can work.

 3. Clone this copy to your local disk (you need the `git` program to do
    this)::

        $ git clone git@github.com:YourLogin/scikit-learn.git

 4. Work on this copy, on your computer, using git to do the version
    control::

        $ git add modified_files
        $ git commit
        $ git push origin master

    and so on.

If your changes are not just trivial fixes, it is better to directly
work in a branch with the name of the feature your are working on. In
this case, replace step 4 by step 5:

  5. Create a branch to host your changes and publish it on your public
     repo::

        $ git checkout -b my-feature
        $ git add modified_files
        $ git commit
        $ git push origin my-feature

When you are ready, and you have pushed your changes on your github repo, go
the web page of the repo, and click on 'Pull request' to send us a pull
request. This will send an email to the commiters, but might also send an
email to the mailing list in order to get more visibility.

It is recommented to check that your contribution complies with the following
rules before submitting a pull request:

    * Follow the `coding-guidelines`_ (see below).

    * When applicable, use the Validation tools and other code in the
      ``sklearn.utils`` submodule.  A list of utility routines available
      for developers can be found in the :ref:`developers-utils` page.

    * All public methods should have informative docstrings with sample
      usage presented as doctests when appropriate.

    * All other tests pass when everything is rebuilt from scrath, under Unix,
      check with (from the toplevel source folder)::

        $ make

    * At least one example script in the ``examples/`` folder. Have a look at
      other examples for reference. Example should demonstrate why this method
      is useful in practice and if possible compare it to other methods
      available in the scikit.

    * At least one paragraph of narrative documentation with links to
      references in the literature (with PDF links when possible) and
      the example.

      The documentation should also include expected time and space
      complexity of the algorithm and scalablity, e.g. "this algorithm can
      scale to a large number of samples > 100000, but does not scale in
      dimensionality: n_features is expected to be lower than 100".

      To build the documentation see `documentation`_ below.

You can also check for common programming errors with the following tools:

    * Code with a good unittest coverage (at least 80%), check with::

        $ pip install nose coverage
        $ nosetests --with-coverage path/to/tests_for_package

    * No pyflakes warnings, check with::

        $ pip install pyflakes
        $ pyflakes path/to/module.py

    * No PEP8 warnings, check with::

        $ pip install pep8
        $ pep8 path/to/module.py

Bonus points for contributions that include a performance analysis with
a benchmark script and profiling output (please report on the mailing
list or on the github wiki).

Also check out the following guide on :ref:`performance-howto` for more
details on profiling and cython optimizations.

.. note::

  The current state of the scikit-learn code base is not compliant with
  all of those guidelines but we expect that enforcing those constraints
  on all new contributions will get the overall code base quality in the
  right direction.

EasyFix Issues
--------------

The best way to get your feet wet is
to pick up an issue from the `issue tracker
<https://github.com/scikit-learn/scikit-learn/issues?labels=EasyFix>`_
that are labeled as EasyFix. This means that the knowledge needed to solve
the issue is low, but still you are helping the project and letting more
experienced developers concentrate on other issues.

.. _contribute_documentation:

Documentation
-------------

We are glad to accept any sort of documentation: function docstrings,
rst docs (like this one), tutorials, etc. Rst docs live in the source
code repository, under directory doc/.

You can edit them using any text editor and generate the html docs by
typing from the doc/ directory ``make html`` (or ``make html-noplot``,
see README in that directory for more info). That should create a
directory _build/html/ with html files that are viewable in a web
browser.

For building the documentation, you will need `sphinx
<http://sphinx.pocoo.org/>`_ and `matplotlib
<http://matplotlib.sourceforge.net/>`_.

When you are writing documentation, it is important to keep a good
compromise between mathematical and algorithmic details, and giving
intuitions to the reader on what the algorithm does. It is best to always
start with a small paragraph with a hand waiving explanation of what the
method does to the data and a figure (coming from an example) ilustrating
it.

.. warning:: **Sphinx version**

   While we do our best to have the documentation build under as many
   version of Sphinx as possible, the different versions tend to behave
   slightly differently. To get the best results, you should use version
   1.0.

Developers web site
-------------------

More information can be found at the `developer's wiki
<https://github.com/scikit-learn/scikit-learn/wiki>`_.


Other ways to contribute
========================

Code is not the only way to contribute to this project. For instance,
documentation is also a very important part of the project and ofter
doesn't get as much attention as it deserves. If you find a typo in
the documentation, or have made improvements, don't hesitate to send
an email to the mailing list or a github pull request. Full
documentation can be found under directory doc/.

It also helps us if you spread the word: reference it from your blog,
articles, link to us from your website, or simply by saying "I use
it":

.. raw:: html

   <script type="text/javascript" src="http://www.ohloh.net/p/480792/widgets/project_users.js?style=rainbow"></script>


.. _coding-guidelines:

Coding guidelines
=================

The following are some guidelines on how new code should be written. Of
course, there are special cases and there will be exceptions to these
rules. However, following these rules when submitting new code makes
the review easier so new code can be integrated in less time.

Uniformly formatted code makes it easier to share code ownership. The
scikit learn tries to follow closely the official Python guidelines
detailed in `PEP8 <http://www.python.org/dev/peps/pep-0008/>`_ that
details how code should be formatted, and indented. Please read it and
follow it.

In addition, we add the following guidelines:

    * Use underscores to separate words in non class names: ``n_samples``
      rather than ``nsamples``.

    * Avoid multiple statements on one line. Prefer a line return after
      a control flow statement (``if``/``for``).

    * Use relative imports for references inside scikit-learn.

    * **Please don't use `import *` in any case**. It is considered harmful
      by the `official Python recommendations
      <http://docs.python.org/howto/doanddont.html#from-module-import>`_.
      It makes the code harder to read as the origin of symbols is no
      longer explicitly referenced, but most important, it prevents
      using a static analysis tool like `pyflakes
      <http://www.divmod.org/trac/wiki/DivmodPyflakes>`_ to automatically
      find bugs in scikit.

    * Use the `numpy docstring standard
      <https://github.com/numpy/numpy/blob/master/doc/HOWTO_DOCUMENT.rst.txt>`_
      in all your docstrings.


A good example of code that we like can be found `here
<https://svn.enthought.com/enthought/browser/sandbox/docs/coding_standard.py>`_.

Input validation
----------------

.. currentmodule:: sklearn.utils

The module :mod:`sklearn.utils` contains various functions for doing input
validation/conversion. Sometimes, ``np.asarray`` suffices for validation;
do `not` use ``np.asanyarray`` or ``np.atleast_2d``, since those let NumPy's
``np.matrix`` through, which has a different API
(e.g., ``*`` means dot product on ``np.matrix``,
but Hadamard product on ``np.ndarray``).

In other cases, be sure to call :func:`safe_asarray`, :func:`atleast2d_or_csr`,
:func:`as_float_array` or :func:`array2d` on any array-like argument passed to a
scikit-learn API function. The exact function to use depends mainly on whether
``scipy.sparse`` matrices must be accepted.

For more information, refer to the :ref:`developers-utils` page.

Random Numbers
--------------

If your code depends on a random number generator, do not use
``numpy.random.random()`` or similar routines.  To ensure
repeatability in error checking, the routine should accept a keyword
``random_state`` and use this to construct a
``numpy.random.RandomState`` object.
See :func:`sklearn.utils.check_random_state` in :ref:`developers-utils`.

Here's a simple example of code using some of the above guidelines:

::

    from sklearn.utils import array2d, check_random_state

    def choose_random_sample(X, random_state=0):
      	"""
	Choose a random point from X

	Parameters
	----------
	X : array-like, shape = (n_samples, n_features)
	    array representing the data
        random_state : RandomState or an int seed (0 by default)
            A random number generator instance to define the state of the
            random permutations generator.

	Returns
	-------
	x : numpy array, shape = (n_features,)
	    A random point selected from X
	"""
    	X = array2d(X)
        random_state = check_random_state(random_state)
	i = random_state.randint(X.shape[0])
	return X[i]


APIs of scikit-learn objects
============================

To have a uniform API, we try to have a common basic API for all the
objects. In addition, to avoid the proliferation of framework code, we
try to adopt simple conventions and limit to a minimum the number of
methods an object has to implement.

Different objects
-----------------

The main objects of the scikit learn are (one class can implement
multiple interfaces):

:Estimator:

    The base object, implements::

      estimator = obj.fit(data)

:Predictor:

    For supervised learning, or some unsupervised problems, implements::

      prediction = obj.predict(data)

:Transformer:

    For filtering or modifying the data, in a supervised or unsupervised
    way, implements::

      new_data = obj.transform(data)

    When fitting and transforming can be performed much more efficiently
    together than separately, implements::

      new_data = obj.fit_transform(data)

:Model:

    A model that can give a goodness of fit or a likelihood of unseen
    data, implements (higher is better)::

      score = obj.score(data)

Estimators
----------

The API has one predominant object: the estimator. A estimator is an
object that fits a model based on some training data and is capable of
inferring some properties on new data. It can be for instance a
classifier or a regressor. All estimators implement the fit method::

    estimator.fit(X, y)

All built-in estimators also have a ``set_params`` method, which sets
data-independent parameters (overriding previous parameter values passed
to ``__init__``). This method is not required for an object to be an
estimator.

All estimators should inherit from ``sklearn.base.BaseEstimator``.

Instantiation
^^^^^^^^^^^^^

This concerns the object creation. The object's ``__init__`` method might
accept as arguments constants that determine the estimator behavior
(like the C constant in SVMs).

It should not, however, take the actual training data as argument, as
this is left to the ``fit()`` method::

    clf2 = SVC(C=2.3)
    clf3 = SVC([[1, 2], [2, 3]], [-1, 1]) # WRONG!


The arguments that go in the ``__init__`` should all be keyword arguments
with a default value. In other words, a user should be able to instanciate
an estimator without passing to it any arguments.

The arguments in given at instanciation of an estimator should all
correspond to hyper parameters describing the model or the optimisation
problem that estimator tries to solve.

In addition, **every keyword argument given to the ``__init__`` should
correspond to an attribute on the instance**. The scikit relies on this
to find what are the relevent attributes to set on an estimator when
doing model selection.

To summarize, a `__init__` should look like::

    def __init__(self, param1=1, param2=2):
        self.param1 = param1
        self.param2 = param2

There should be no logic, and the parameters should not be changed.
The corresponding logic should be put when the parameters are used. The
following is wrong::

    def __init__(self, param1=1, param2=2, param3=3):
        # WRONG: parameters should not be modified
        if param1 > 1:
            param2 += 1
        self.param1 = param1
        # WRONG: the object's attributes should have exactly the name of
        # the argument in the constructor
        self.param3 = param2

Scikit-Learn relies on this mechanism to introspect object to set
their parameters by cross-validation.

Fitting
^^^^^^^

The next thing you'll probably want to do is to estimate some
parameters in the model. This is implemented in the .fit() method.

The fit method takes as argument the training data, which can be one
array in the case of unsupervised learning, or two arrays in the case
of supervised learning.

Note that the model is fitted using X and y but the object holds no
reference to X, y. There are however some exceptions to this, as in
the case of precomputed kernels where you need to store access these
data in the predict method.

============= ======================================================
Parameters
============= ======================================================
X             array-like, with shape = [N, D], where N is the number
              of samples and D is the number of features.

y             array, with shape = [N], where N is the number of
              samples.

kwargs        optional data dependent parameters.
============= ======================================================

``X.shape[0]`` should be the same as ``y.shape[0]``. If this requisite
is not met, an exception of type ``ValueError`` should be raised.

``y`` might be ignored in the case of unsupervised learning. However to
make it possible to use the estimator as part of a pipeline that can
mix both supervised and unsupervised transformers even unsupervised
estimators are kindly ask to accept a ``y=None`` keyword argument in
the second position that is just ignored by the estimator.

The method should return the object (``self``). This pattern is useful
to be able to implement quick one liners in an ipython session such as::

  y_predicted = SVC(C=100).fit(X_train, y_train).predict(X_test)

Depending on the nature of the algorithm ``fit`` can sometimes also
accept additional keywords arguments. However any parameter that can
have a value assigned prior having access to the data should be an
``__init__`` keyword argument. **fit parameters should be restricted
to directly data dependent variables**. For instance a Gram matrix or
an affinity matrix which are precomputed from the data matrix ``X`` are
data dependent. A tolerance stopping criterion ``tol`` is not directly
data dependent (although the optimal value according to some scoring
function probably is).

Any attribute that ends with ``_`` is expected to be overridden when
you call ``fit`` a second time without taking any previous value into
account: **fit should be idempotent**.

Optional Arguments
^^^^^^^^^^^^^^^^^^

In iterative algorithms, number of iterations should be specified by
an int called ``n_iter``.

Unresolved API issues
----------------------

Some things are must still be decided:

    * what should happen when predict is called before than fit() ?
    * which exception should be raised when arrays' shape do not match
      in fit() ?

Working notes
---------------

For unresolved issues, TODOs, remarks on ongoing work, developers are
adviced to maintain notes on the github wiki:
https://github.com/scikit-learn/scikit-learn/wiki

Specific models
-----------------

In linear models, coefficients are stored in an array called ``coef_``,
and independent term is stored in ``intercept_``.
