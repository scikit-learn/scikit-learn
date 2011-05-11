============
Contributing
============

This project is a community effort, and everyone is welcomed to
contribute.

The project is hosted on http://github.com/scikit-learn/scikit-learn

Submitting a bug report
=======================

In case you experience difficulties using the package, do not hesitate
to submit a ticket to the
`Bug Tracker <http://github.com/scikit-learn/scikit-learn/issues>`_.

You are also welcomed to post there feature requests and patches.

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

.. note:

  Before to starting to work on a non trivial new feature, it highly advised
  to discuss it on the developer mailing list.

  https://lists.sourceforge.net/lists/listinfo/scikit-learn-general

  The goal is to avoid duplicated work (this has occurred several times in the
  past).


How to contribute
-----------------

The prefered way to contribute to `scikit-learn` is to fork the main
repository on
`github <http://github.com/scikit-learn/scikit-learn/>`__:

 1. `Create an account <https://github.com/signup/free>`_ on
    github if you don't have one already.

 2. Fork the `scikit-learn repo
    <http://github.com/scikit-learn/scikit-learn>`__: click on the 'Fork'
    button, at the top, center of the page. This creates a copy of
    the code on the github server where you can work.

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
        $ gid add modified_files
        $ git commit
        $ git push origin my-feature

When you are ready, and you have pushed your changes on your github repo,
go the web page of the repo, and click on 'Pull request' to send us a
pull request. Send us a mail with your pull request, and we can look at
your changes, and integrate them.

**Before asking for a pull or a review**, please check that your contribution
complies with the following rules:

    * Follow the `coding-guidelines`_ (see below).

    * All public methods should have informative docstrings with sample
      usage presented as doctests when appropriate.

    * Code with a good unittest coverage (at least 80%), check with::

        $ pip install nose coverage
        $ nosetests --with-coverage path/to/tests_for_package

    * All other tests pass when everything is rebuilt from scrath, under Unix,
      check with (from the toplevel source folder)::

        $ make

    * No pyflakes warnings, check with::

        $ pip install pyflakes
        $ pyflakes path/to/module.py

    * No PEP8 warnings, check with::

        $ pip install pep8
        $ pep8 path/to/module.py

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

    * Use relative imports for references inside scikits.learn.

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


APIs of scikit learn objects
=============================

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


Instantiation
^^^^^^^^^^^^^^

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
problem that estimator tries to solve. They should however not be
parameters of the estimation routine: these are passed directly to the
``fit`` method.

In addition, **every keyword argument given to the ``__init__`` should
correspond to an attribute on the instance**. The scikit relies on this
to find what are the relevent attributes to set on an estimator when
doing model selection.

All estimators should inherit from ``scikit.learn.base.BaseEstimator``.


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

Y             array, with shape = [N], where N is the number of 
              samples.

args, kwargs  Parameters can also be set in the fit method.
============= ======================================================

X.shape[0] should be the same as Y.shape[0]. If this requisite is not
met, an exception should be raised.

Y might be dropped in the case of unsupervised learning.

The method should return the object (``self``).


Python tuples
^^^^^^^^^^^^^^

In addition to numpy arrays, all methods should be able to accept
Python tuples as arguments. In practice, this means you should call
``numpy.asanyarray`` at the beginning at each public method that accepts
arrays.


Optional Arguments
^^^^^^^^^^^^^^^^^^^

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
