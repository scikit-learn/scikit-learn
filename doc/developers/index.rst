===============
Contributing
===============

This project is a community effort, and everyone is welcomed to
contribute.


Submitting a bug report 
=========================

In case you experience difficulties using the package, do not hesitate
to submit a ticket to the
`Bug Tracker <http://sourceforge.net/apps/trac/scikit-learn/report/1>`_.

You are also welcomed to post there feature requests and patches.

.. _git_repo:

Retrieving the latest code
==========================

You can check the latest sources with the command::

    git clone git://github.com/scikit-learn/scikit-learn.git

or if you have write privileges::

    git clone git@github.com:scikit-learn/scikit-learn.git

You can also check out the sources online in the web page
http://github.com/scikit-learn/scikit-learn 

If you run the development version, it is cumbersome to re-install the
package each time you update the sources. It is thus preferred that
you add the scikit-directory to your PYTHONPATH and build the
extension in place::

    python setup.py build_ext --inplace


Contributing code
===========================

How to contribute
-------------------

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

When you are ready, and you have pushed your changes on your github repo,
go the web page of the repo, and click on 'Pull request' to send us a
pull request. Send us a mail with your pull request, and we can look at
your changes, and integrate them.

**Before asking for a pull or a review**, be sure to read the 
`coding-guidelines`_ (below).

Also, make sure that your code is tested, and that all the tests for the
scikit pass.

EasyFix Issues
---------------

The best way to get your feet wet is to pick up an issue from the
`issue tracker
<https://sourceforge.net/apps/trac/scikit-learn/report>`_ that are
labeled as EasyFix. This means that the knowledge needed to solve the
issue is low, but still you are helping the project and letting more
experienced developers concentrate on other issues.


Roadmap
-------

`Here <http://sourceforge.net/apps/trac/scikit-learn/roadmap>`_ you
will find a detailed roadmap, with a description on what's planned to
be implemented in the following releases.


Documentation
----------------------

We are glad to accept any sort of documentation: function docstrings,
rst docs (like this one), tutorials, etc. Rst docs live in the source
code repository, under directory doc/.

You can edit them using any text editor and generate the html docs by
typing from the doc/ directory ``make html`` (or ``make html-noplot``,
see README in that directory for more info). That should create a
directory _build/html/ with html files that are viewable in a web
browser.


Developers web site
----------------------

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
===================

The following are some guidelines on how new code should be
written. Of course, there are special cases and there will be
exceptions to these rules. However, following these rules when
submitting new code makes the review easier so new code can be
integrated in less time.

Uniformly formated code makes it easier to share code ownership. The
scikit learn tries to follow closely the officiel Python guidelines
detailed in `PEP8 <http://www.python.org/dev/peps/pep-0008/>`_ that
details how code should be formatted, and indented. Please read it and
follow it.

In addition, we add the following guidelines:

    * Use underscores to separate words in non class names: `n_samples`
      rather than `nsamples`.

    * Avoid multiple statements on one line. Prefer a line return after
      a control flow statement (`if`/`for`).

    * Use relative imports for references inside scikits.learn.

    * **Please don't use `import *` in any case**. It is considered harmful 
      by the `official Python recommandations
      <http://docs.python.org/howto/doanddont.html#from-module-import>`_.
      It makes the code harder to read as the origin of symbols is no 
      longer explicitely referenced, but most important, it prevents
      using a static analysis tool like `pyflakes
      <http://www.divmod.org/trac/wiki/DivmodPyflakes>`_ to automatically
      find bugs in the scikit.

A good example of code that we like can be found `here
<https://svn.enthought.com/enthought/browser/sandbox/docs/coding_standard.py>`_.

APIs of scikit learn objects
=============================

To have a uniform API, we try to have a common basic API for all the
objects. In addition, to avoid the proliferation of framework code, we
try to adopt simple conventions and limit to a minimum the number of
methods an object has to implement.

Different objects
-------------------

The main objects of the scikit learn are (one class can implement
multiple interfaces):

:Estimator:

    The base object, implements::

	estimator = obj.fit(data)

:Predictor:

    For suppervised learning, or some unsupervised problems, implements::

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
--------------

The API has one predominant object: the estimator. A estimator is an
object that fits a model based on some training data and is capable of
inferring some properties on new data. It can be for instance a
classifier or a regressor. All estimators implement the fit method::

    estimator.fit(X, y)


Instantiation
^^^^^^^^^^^^^^

This concerns the object creation. The object's __init__ method might
accept as arguments constants that determine the estimator behavior
(like the C constant in SVMs).

It should not, however, take the actual training data as argument, as
this is left to the ``fit()`` method::

    clf2 = SVC(C=2.3)
    clf3 = SVC([[1, 2], [2, 3]], [-1, 1]) # WRONG!


The arguments that go in the `__init__` should all be keyword arguments
with a default value. In other words, a user should be able to instanciate
an estimator without passing to it any arguments.

The arguments in given at instanciation of an estimator should all
correspond to hyper parameters describing the model or the optimisation
problem that estimator tries to solve. They should however not be
parameters of the estimation routine: these are passed directly to the
`fit` method. 

In addition, **every keyword argument given to the `__init__` should
correspond to an attribute on the instance**. The scikit relies on this
to find what are the relevent attributes to set on an estimator when
doing model selection.

All estimators should inherit from `scikit.learn.base.BaseEstimator`

Fitting
^^^^^^^^^^^^^^

The next thing you'll probably want to do is to estimate some
parameters in the model. This is implemented in the .fit() method.

The fit method takes as argument the training data, which can be one
array in the case of unsupervised learning, or two arrays in the case
of supervised learning.

Note that the model is fitted using X and y but the object holds no
reference to X, y. There are however some exceptions to this, as in
the case of precomputed kernels where you need to store access these
data in the predict method.

  Parameters

    * X : array-like, with shape = [N, D], where N is the number of
      samples and D is the number of features.
    * Y : array, with shape = [N], where N is the number of samples.

    * args, kwargs. Parameters can also be set in the fit method.

X.shape[0] should be the same as Y.shape[0]. If this requisite is not
met, an exception should be raised.

Y might be dropped in the case of unsupervised learning.

The method should return the object (self).


Python tuples
^^^^^^^^^^^^^^

In addition to numpy arrays, all methods should be able to accept
python tuples as arguments. In practice, this means you should call
numpy.asanyarray at the beginning at each public method that accepts
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


Specific models
-----------------

In linear models, coefficients are stored in an array called ``coef_``,
and independent term is stored in ``intercept_``.
