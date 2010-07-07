===========
Development
===========

This project is a community effort, and everyone is welcomed to
contribute.

Bug Tracker
===========

In case you experience difficulties using the package, do not hesitate
to submit a ticket to the
`Bug Tracker <http://sourceforge.net/apps/trac/scikit-learn/report/1>`_.

You are also welcomed to post there feature requests and patches.

Code
====


Git repo
--------

You can check the latest sources with the command::

    git clone git://scikit-learn.git.sourceforge.net/gitroot/scikit-learn/scikit-learn

or if you have write privileges::

    git clone ssh://USERNAME@scikit-learn.git.sourceforge.net/gitroot/scikit-learn/scikit-learn

If you have contributed some code and would like to have write
privileges in subversion repository, please contact me (Fabian
Pedregosa <fabian.pedregosa@inria.fr>) and I'll give you write
privileges.

If you run the development version, it is cumbersome to re-install the
package each time you update the sources. It is thus preferred that
you add the scikit-directory to your PYTHONPATH and build the
extension in place::

    python setup.py build_ext --inplace


Patches
-------
Patches are the prefered way to contribute to a project if you do not
have write privileges.

Let's suppose that you have the latest sources for subversion and that
you just made some modifications that you'd like to share with the
world. You might proceed as:

1. Create a patch file. The command::

    git format-patch origin

will create a series of patch files with the changes you made with
the code base. 

2. Send that file to the mailing list or attach it to an
issue in the issue tracker and some devs will push that patch to the
main repository.

3. Wait for a reply. You should soon receive a reply on whether your
patch was committed.


EasyFix Issues
^^^^^^^^^^^^^^

The best way to get your feet wet is to pick up an issue from the
`issue tracker
<https://sourceforge.net/apps/trac/scikit-learn/report>`_ that are
labeled as EasyFix. This means that the knowledge needed to solve the
issue is low, but still you are helping the project and letting more
experienced developers concentrate on other issues.


Roadmap
-------

`Here <http://sourceforge.net/apps/trac/scikit-learn/roadmap`_ you
will find a detailed roadmap, with a description on what's planned to
be implemented in the following releases.

.. _packaging:

Packaging
^^^^^^^^^

You can also help making binary distributions for windows, OsX or packages for some
distribution.

Developers web site
===================
More information can be found at the developer's web site:
http://sourceforge.net/apps/trac/scikit-learn/wiki , which contains a
wiki, an issue tracker, and a Roadmap

Documentation
=============

I am glad to accept any sort of documentation: function docstrings,
rst docs (like this one), tutorials, etc. Rst docs live in the source
code repository, under directory doc/.

You can edit them using any text editor and generate the html docs by
typing ``make html`` from the doc/ directory. That should create a
directory _build/html/ with html files that are viewable in a web
browser.


API guidelines
==============

The following are some guidelines on how new code should be
written. Of course, there are special cases and there will be
exceptions to these rules. However, following these rules when
submitting new code makes the review easier so new code can be
integrated in less time.


Estimators
----------

The API has one predominant object: the estimator. A estimator is an
object that fits a model based on some training data and is capable of
inferring some properties on new data. It can be for instance a
classifier or a regressor.


Instantiation
^^^^^^^^^^^^^

This concerns the object creation. The object's __init__ method might
accept as arguments constants that determine the estimator behavior
(like the C constant in SVMs).

It should not, however, take the actual training data as argument, as
this is leaved to the ``fit()`` method::

    clf1 = SVM(impl='c_svm')
    clf2 = SVM(C=2.3)
    clf3 = SVM([[1, 2], [2, 3]], [-1, 1]) # WRONG!


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

  Parameters

    * X : array-like, with shape = [N, D], where N is the number of
      samples and D is the number of features.
    * Y : array, with shape = [N], where N is the number of samples.

X.shape[0] should be the same as Y.shape[0]. If this requisite is not
met, an exception should be raised.

Y might be dropped in the case of unsupervised learning.

The method should return the object (self).


Python tuples
^^^^^^^^^^^^^

In addition to numpy arrays, all methods should be able to accept
python tuples as arguments. In practice, this means you should call
numpy.asanyarray at the beginning at each public method that accepts
arrays.


Optional Arguments
^^^^^^^^^^^^^^^^^^

In iterative algorithms, number of iterations should be specified by
an int called ``n_iter``.


TODO
----
Some things are must still be decided:

    * what should happen when predict is called before than fit() ?
    * which exception should be raised when arrays' shape do not match
      in fit() ?


Specific models
---------------

In linear models, coefficients are stored in an array called ``coef_``,
and independent term is stored in ``intercept_``.
