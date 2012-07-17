.. _contributing:

============
Contributing
============

This project is a community effort, and everyone is welcome to
contribute.

The project is hosted on http://github.com/scikit-learn/scikit-learn


Submitting a bug report
=======================

In case you experience issues using this package, do not hesitate to submit a
ticket to the
`Bug Tracker <http://github.com/scikit-learn/scikit-learn/issues>`_. You are
also welcome to post feature requests or links to pull requests.


.. _git_repo:

Retrieving the latest code
==========================

We use `Git <http://git-scm.com/>`_ for version control and
`GitHub <http://github.com/>`_ for hosting our main repository.

You can check out the latest sources with the command::

    git clone git://github.com/scikit-learn/scikit-learn.git

or if you have write privileges::

    git clone git@github.com:scikit-learn/scikit-learn.git

If you run the development version, it is cumbersome to reinstall the
package each time you update the sources. It is thus preferred that
you add the scikit-learn directory to your ``PYTHONPATH`` and build the
extension in place::

    python setup.py build_ext --inplace

On Unix-like systems, you can simply type ``make`` in the top-level folder to
build in-place and launch all the tests. Have a look at the ``Makefile`` for
additional utilities.


Contributing code
=================

.. note::

  To avoid duplicating work, it is highly advised that you contact the
  developers on the mailing list before starting work on a non-trivial feature.

  https://lists.sourceforge.net/lists/listinfo/scikit-learn-general

How to contribute
-----------------

The preferred way to contribute to scikit-learn is to fork the `main
repository <http://github.com/scikit-learn/scikit-learn/>`__ on GitHub:

 1. `Create an account <https://github.com/signup/free>`_ on
    GitHub if you do not already have one.

 2. Fork the `project repository
    <http://github.com/scikit-learn/scikit-learn>`__: click on the 'Fork'
    button near the top of the page. This creates a copy of the code under your
    account on the GitHub server.

 3. Clone this copy to your local disk::

        $ git clone git@github.com:YourLogin/scikit-learn.git

 4. Create a branch to hold your changes::

        $ git checkout -b my-feature

    and start making changes. Never work in the ``master`` branch!

 5. Work on this copy, on your computer, using Git to do the version
    control. When you're done editing, do::

        $ git add modified_files
        $ git commit

    to record your changes in Git, then push them to GitHub with::

        $ git push -u origin my-feature

Finally, go to the web page of the your fork of the scikit-learn repo,
and click 'Pull request' to send your changes to the maintainers for review.
request. This will send an email to the committers, but might also send an
email to the mailing list in order to get more visibility.

.. note::

  In the above setup, your ``origin`` remote repository points to
  YourLogin/scikit-learn.git. If you wish to `fetch/merge` from the main
  repository instead of your `forked` one, you will need to add another remote
  to use instead of ``origin``. If we choose the name ``upstream`` for it, the
  command will be::

        $ git remote add upstream https://github.com/scikit-learn/scikit-learn.git

(If any of the above seems like magic to you, then look up the
`Git documentation <http://git-scm.com/documentation>`_ on the web.)

It is recommended to check that your contribution complies with the following
rules before submitting a pull request:

    * Follow the `coding-guidelines`_ (see below).

    * When applicable, use the Validation tools and other code in the
      ``sklearn.utils`` submodule.  A list of utility routines available
      for developers can be found in the :ref:`developers-utils` page.

    * All public methods should have informative docstrings with sample
      usage presented as doctests when appropriate.

    * All other tests pass when everything is rebuilt from scratch. On
      Unix-like systems, check with (from the toplevel source folder)::

        $ make

    * When adding additional functionality, provide at least one example script
      in the ``examples/`` folder. Have a look at other examples for reference.
      Examples should demonstrate why the new functionality is useful in
      practice and, if possible, compare it to other methods available in
      scikit-learn.

    * At least one paragraph of narrative documentation with links to
      references in the literature (with PDF links when possible) and
      the example.

      The documentation should also include expected time and space
      complexity of the algorithm and scalability, e.g. "this algorithm can
      scale to a large number of samples > 100000, but does not scale in
      dimensionality: n_features is expected to be lower than 100".

      To build the documentation, see the `documentation`_ section below.

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

    * AutoPEP8 can help you fix some of the easy redundant errors::

        $ pip install autopep8
        $ autopep8 path/to/pep8.py

Bonus points for contributions that include a performance analysis with
a benchmark script and profiling output (please report on the mailing
list or on the GitHub wiki).

Also check out the :ref:`performance-howto` guide for more details on profiling
and Cython optimizations.

.. note::

  The current state of the scikit-learn code base is not compliant with
  all of those guidelines, but we expect that enforcing those constraints
  on all new contributions will get the overall code base quality in the
  right direction.

EasyFix Issues
--------------

A great way to start contributing to scikit-learn is to pick an item from the
list of `EasyFix issues
<https://github.com/scikit-learn/scikit-learn/issues?labels=EasyFix>`_
in the issue tracker.  Resolving these issues allow you to start contributing
to the project without much prior knowledge. Your assistance in this area will
be greatly appreciated by the more experienced developers as it helps free up
their time to concentrate on other issues.

.. _contribute_documentation:

Documentation
-------------

We are glad to accept any sort of documentation: function docstrings,
reStructuredText documents (like this one), tutorials, etc. reStructuredText
documents live in the source code repository under the doc/ directory.

You can edit the documentation using any text editor, and then generate the
HTML output by typing ``make html`` from the doc/ directory. Alternatively,
``make html-noplot`` can be used to quickly generate the documentation without
the example gallery. The resulting HTML files will be placed in _build/html/
and are viewable in a web browser. See the README file in the doc/ directory
for more information.

For building the documentation, you will need `sphinx
<http://sphinx.pocoo.org/>`_ and `matplotlib
<http://matplotlib.sourceforge.net/>`_.

When you are writing documentation, it is important to keep a good
compromise between mathematical and algorithmic details, and give
intuition to the reader on what the algorithm does. It is best to always
start with a small paragraph with a hand-waiving explanation of what the
method does to the data and a figure (coming from an example) illustrating
it.

.. warning:: **Sphinx version**

   While we do our best to have the documentation build under as many
   version of Sphinx as possible, the different versions tend to behave
   slightly differently. To get the best results, you should use version
   1.0.

Developers web site
-------------------

More information can be found on the `developer's wiki
<https://github.com/scikit-learn/scikit-learn/wiki>`_.


Other ways to contribute
========================

Code is not the only way to contribute to scikit-learn. For instance,
documentation is also a very important part of the project and often
doesn't get as much attention as it deserves. If you find a typo in
the documentation, or have made improvements, do not hesitate to send
an email to the mailing list or submit a GitHub pull request. Full
documentation can be found under the doc/ directory.

It also helps us if you spread the word: reference the project from your blog
and articles, link to it from your website, or simply say "I use it":

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
scikit-learn project tries to closely follow the official Python guidelines
detailed in `PEP8 <http://www.python.org/dev/peps/pep-0008/>`_ that
detail how code should be formatted and indented. Please read it and
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
      find bugs in scikit-learn.

    * Use the `numpy docstring standard
      <https://github.com/numpy/numpy/blob/master/doc/HOWTO_DOCUMENT.rst.txt>`_
      in all your docstrings.


A good example of code that we like can be found `here
<https://svn.enthought.com/enthought/browser/sandbox/docs/coding_standard.py>`_.

Input validation
----------------

.. currentmodule:: sklearn.utils

The module :mod:`sklearn.utils` contains various functions for doing input
validation and conversion. Sometimes, ``np.asarray`` suffices for validation;
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

Here's a simple example of code using some of the above guidelines::

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
methods an object must implement.

Different objects
-----------------

The main objects in scikit-learn are (one class can implement
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
inferring some properties on new data. It can be, for instance, a
classifier or a regressor. All estimators implement the fit method::

    estimator.fit(X, y)

All built-in estimators also have a ``set_params`` method, which sets
data-independent parameters (overriding previous parameter values passed
to ``__init__``). This method is not required for an object to be an
estimator.

All estimators should inherit from ``sklearn.base.BaseEstimator``.

Instantiation
^^^^^^^^^^^^^

This concerns the creation of an object. The object's ``__init__`` method
might accept constants as arguments that determine the estimator's behavior
(like the C constant in SVMs). It should not, however, take the actual training
data as an argument, as this is left to the ``fit()`` method::

    clf2 = SVC(C=2.3)
    clf3 = SVC([[1, 2], [2, 3]], [-1, 1]) # WRONG!


The arguments accepted by ``__init__`` should all be keyword arguments
with a default value. In other words, a user should be able to instantiate
an estimator without passing any arguments to it. The arguments should all
correspond to hyperparameters describing the model or the optimisation
problem the estimator tries to solve.

In addition, **every keyword argument accepted by ``__init__`` should
correspond to an attribute on the instance**. Scikit-learn relies on this to
find the relevant attributes to set on an estimator when doing model selection.

To summarize, a `__init__` should look like::

    def __init__(self, param1=1, param2=2):
        self.param1 = param1
        self.param2 = param2

There should be no logic, and the parameters should not be changed.
The corresponding logic should be put where the parameters are used. The
following is wrong::

    def __init__(self, param1=1, param2=2, param3=3):
        # WRONG: parameters should not be modified
        if param1 > 1:
            param2 += 1
        self.param1 = param1
        # WRONG: the object's attributes should have exactly the name of
        # the argument in the constructor
        self.param3 = param2

Scikit-learn relies on this mechanism to introspect objects to set
their parameters by cross-validation.

Fitting
^^^^^^^

The next thing you will probably want to do is to estimate some
parameters in the model. This is implemented in the ``fit()`` method.

The ``fit()`` method takes the training data as arguments, which can be one
array in the case of unsupervised learning, or two arrays in the case
of supervised learning.

Note that the model is fitted using X and y, but the object holds no
reference to X and y. There are, however, some exceptions to this, as in
the case of precomputed kernels where this data must be stored for use by
the predict method.

============= ======================================================
Parameters
============= ======================================================
X             array-like, with shape = [N, D], where N is the number
              of samples and D is the number of features.

y             array, with shape = [N], where N is the number of
              samples.

kwargs        optional data-dependent parameters.
============= ======================================================

``X.shape[0]`` should be the same as ``y.shape[0]``. If this requisite
is not met, an exception of type ``ValueError`` should be raised.

``y`` might be ignored in the case of unsupervised learning. However, to
make it possible to use the estimator as part of a pipeline that can
mix both supervised and unsupervised transformers, even unsupervised
estimators are kindly asked to accept a ``y=None`` keyword argument in
the second position that is just ignored by the estimator.

The method should return the object (``self``). This pattern is useful
to be able to implement quick one liners in an IPython session such as::

  y_predicted = SVC(C=100).fit(X_train, y_train).predict(X_test)

Depending on the nature of the algorithm, ``fit`` can sometimes also
accept additional keywords arguments. However, any parameter that can
have a value assigned prior to having access to the data should be an
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

In iterative algorithms, the number of iterations should be specified by
an integer called ``n_iter``.

Unresolved API issues
----------------------

Some things are must still be decided:

    * what should happen when predict is called before ``fit()`` ?
    * which exception should be raised when the shape of arrays do not match
      in ``fit()`` ?

Working notes
-------------

For unresolved issues, TODOs, and remarks on ongoing work, developers are
advised to maintain notes on the `GitHub wiki
<https://github.com/scikit-learn/scikit-learn/wiki>`__.

Specific models
---------------

In linear models, coefficients are stored in an array called ``coef_``,
and the independent term is stored in ``intercept_``.
