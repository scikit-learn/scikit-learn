.. _develop:

==================================
Developing scikit-learn estimators
==================================

Whether you are proposing an estimator for inclusion in scikit-learn,
developing a separate package compatible with scikit-learn, or
implementing custom components for your own projects, this chapter
details how to develop objects that safely interact with scikit-learn
pipelines and model selection tools.

This section details the public API you should use and implement for a scikit-learn
compatible estimator. Inside scikit-learn itself, we experiment and use some private
tools and our goal is always to make them public once they are stable enough, so that
you can also use them in your own projects.

.. currentmodule:: sklearn

.. _api_overview:

APIs of scikit-learn objects
============================

There are two major types of estimators. You can think of the first group as simple
estimators, which consists most estimators, such as
:class:`~sklearn.linear_model.LogisticRegression` or
:class:`~sklearn.ensemble.RandomForestClassifier`. And the second group are
meta-estimators, which are estimators that wrap other estimators.
:class:`~sklearn.pipeline.Pipeline` and :class:`~sklearn.model_selection.GridSearchCV`
are two examples of meta-estimators.

Here we start with a few vocabulary, and then we illustrate how you can implement
your own estimators.

Elements of the scikit-learn API are described more definitively in the
:ref:`glossary`.

Different objects
-----------------

The main objects in scikit-learn are (one class can implement multiple interfaces):

:Estimator:

    The base object, implements a ``fit`` method to learn from data, either::

      estimator = estimator.fit(data, targets)

    or::

      estimator = estimator.fit(data)

:Predictor:

    For supervised learning, or some unsupervised problems, implements::

      prediction = predictor.predict(data)

    Classification algorithms usually also offer a way to quantify certainty
    of a prediction, either using ``decision_function`` or ``predict_proba``::

      probability = predictor.predict_proba(data)

:Transformer:

    For modifying the data in a supervised or unsupervised way (e.g. by adding, changing,
    or removing columns, but not by adding or removing rows). Implements::

      new_data = transformer.transform(data)

    When fitting and transforming can be performed much more efficiently
    together than separately, implements::

      new_data = transformer.fit_transform(data)

:Model:

    A model that can give a `goodness of fit
    <https://en.wikipedia.org/wiki/Goodness_of_fit>`_ measure or a likelihood of
    unseen data, implements (higher is better)::

      score = model.score(data)

Estimators
----------

The API has one predominant object: the estimator. An estimator is an
object that fits a model based on some training data and is capable of
inferring some properties on new data. It can be, for instance, a
classifier or a regressor. All estimators implement the fit method::

    estimator.fit(X, y)

Out of all the methods that an estimator implements, ``fit`` is usually the one you
want to implement yourself. Other methods such as ``set_params``, ``get_params``, etc.
are implemented in :class:`~sklearn.base.BaseEstimator`, which you should inherit from.
You might need to inherit from more mixins, which we will explain later.

Instantiation
^^^^^^^^^^^^^

This concerns the creation of an object. The object's ``__init__`` method might accept
constants as arguments that determine the estimator's behavior (like the ``alpha``
constant in :class:`~sklearn.linear_model.SGDClassifier`). It should not, however, take
the actual training data as an argument, as this is left to the ``fit()`` method::

    clf2 = SGDClassifier(alpha=2.3)
    clf3 = SGDClassifier([[1, 2], [2, 3]], [-1, 1]) # WRONG!


Ideally, the arguments accepted by ``__init__`` should all be keyword arguments with a
default value. In other words, a user should be able to instantiate an estimator without
passing any arguments to it. In some cases, where there are no sane defaults for an
argument, they can be left without a default value. In scikit-learn itself, we have
very few places, only in some meta-estimators, where the sub-estimator(s) argument is
a required argument.

Most arguments correspond to hyperparameters describing the model or the optimisation
problem the estimator tries to solve. Other parameters might define how the estimator
behaves, e.g. defining the location of a cache to store some data. These initial
arguments (or parameters) are always remembered by the estimator. Also note that they
should not be documented under the "Attributes" section, but rather under the
"Parameters" section for that estimator.

In addition, **every keyword argument accepted by** ``__init__`` **should
correspond to an attribute on the instance**. Scikit-learn relies on this to
find the relevant attributes to set on an estimator when doing model selection.

To summarize, an ``__init__`` should look like::

    def __init__(self, param1=1, param2=2):
        self.param1 = param1
        self.param2 = param2

There should be no logic, not even input validation, and the parameters should not be
changed; which also means ideally they should not be mutable objects such as lists or
dictionaries. If they're mutable, they should be copied before being modified. The
corresponding logic should be put where the parameters are used, typically in ``fit``.
The following is wrong::

    def __init__(self, param1=1, param2=2, param3=3):
        # WRONG: parameters should not be modified
        if param1 > 1:
            param2 += 1
        self.param1 = param1
        # WRONG: the object's attributes should have exactly the name of
        # the argument in the constructor
        self.param3 = param2

The reason for postponing the validation is that if ``__init__`` includes input
validation, then the same validation would have to be performed in ``set_params``, which
is used in algorithms like :class:`~sklearn.model_selection.GridSearchCV`.

Also it is expected that parameters with trailing ``_`` are **not to be set
inside the** ``__init__`` **method**. More details on attributes that are not init
arguments come shortly.

Fitting
^^^^^^^

The next thing you will probably want to do is to estimate some parameters in the model.
This is implemented in the ``fit()`` method, and it's where the training happens.
For instance, this is where you have the computation to learn or estimate coefficients
for a linear model.

The ``fit()`` method takes the training data as arguments, which can be one
array in the case of unsupervised learning, or two arrays in the case
of supervised learning. Other metadata that come with the training data, such as
``sample_weight``, can also be passed to ``fit`` as keyword arguments.

Note that the model is fitted using ``X`` and ``y``, but the object holds no
reference to ``X`` and ``y``. There are, however, some exceptions to this, as in
the case of precomputed kernels where this data must be stored for use by
the predict method.

============= ======================================================
Parameters
============= ======================================================
X             array-like of shape (n_samples, n_features)

y             array-like of shape (n_samples,)

kwargs        optional data-dependent parameters
============= ======================================================

The number of samples, i.e. ``X.shape[0]`` should be the same as ``y.shape[0]``. If this
requirement is not met, an exception of type ``ValueError`` should be raised.

``y`` might be ignored in the case of unsupervised learning. However, to
make it possible to use the estimator as part of a pipeline that can
mix both supervised and unsupervised transformers, even unsupervised
estimators need to accept a ``y=None`` keyword argument in
the second position that is just ignored by the estimator.
For the same reason, ``fit_predict``, ``fit_transform``, ``score``
and ``partial_fit`` methods need to accept a ``y`` argument in
the second place if they are implemented.

The method should return the object (``self``). This pattern is useful
to be able to implement quick one liners in an IPython session such as::

  y_predicted = SGDClassifier(alpha=10).fit(X_train, y_train).predict(X_test)

Depending on the nature of the algorithm, ``fit`` can sometimes also accept additional
keywords arguments. However, any parameter that can have a value assigned prior to
having access to the data should be an ``__init__`` keyword argument. Ideally, **fit
parameters should be restricted to directly data dependent variables**. For instance a
Gram matrix or an affinity matrix which are precomputed from the data matrix ``X`` are
data dependent. A tolerance stopping criterion ``tol`` is not directly data dependent
(although the optimal value according to some scoring function probably is).

When ``fit`` is called, any previous call to ``fit`` should be ignored. In
general, calling ``estimator.fit(X1)`` and then ``estimator.fit(X2)`` should
be the same as only calling ``estimator.fit(X2)``. However, this may not be
true in practice when ``fit`` depends on some random process, see
:term:`random_state`. Another exception to this rule is when the
hyper-parameter ``warm_start`` is set to ``True`` for estimators that
support it. ``warm_start=True`` means that the previous state of the
trainable parameters of the estimator are reused instead of using the
default initialization strategy.

Estimated Attributes
^^^^^^^^^^^^^^^^^^^^

According to scikit-learn conventions, attributes which you'd want to expose to your
users as public attributes and have been estimated or learned from the data must always
have a name ending with trailing underscore, for example the coefficients of some
regression estimator would be stored in a ``coef_`` attribute after ``fit`` has been
called. Similarly, attributes that you learn in the process and you'd like to store yet
not expose to the user, should have a leading underscore, e.g. ``_intermediate_coefs``.
You'd need to document the first group (with a trailing underscore) as "Attributes" and
no need to document the second group (with a leading underscore).

The estimated attributes are expected to be overridden when you call ``fit`` a second
time.

Universal attributes
^^^^^^^^^^^^^^^^^^^^

Estimators that expect tabular input should set a `n_features_in_`
attribute at `fit` time to indicate the number of features that the estimator
expects for subsequent calls to :term:`predict` or :term:`transform`.
See `SLEP010
<https://scikit-learn-enhancement-proposals.readthedocs.io/en/latest/slep010/proposal.html>`__
for details.

Similarly, if estimators are given dataframes such as pandas or polars, they should
set a ``feature_names_in_`` attribute to indicate the features names of the input data,
detailed in `SLEP007
<https://scikit-learn-enhancement-proposals.readthedocs.io/en/latest/slep007/proposal.html>`__.
Using :func:`~sklearn.utils.validation.validate_data` would automatically set these
attributes for you.

.. _rolling_your_own_estimator:

Rolling your own estimator
==========================
If you want to implement a new estimator that is scikit-learn compatible, there are
several internals of scikit-learn that you should be aware of in addition to
the scikit-learn API outlined above. You can check whether your estimator
adheres to the scikit-learn interface and standards by running
:func:`~sklearn.utils.estimator_checks.check_estimator` on an instance. The
:func:`~sklearn.utils.estimator_checks.parametrize_with_checks` pytest
decorator can also be used (see its docstring for details and possible
interactions with `pytest`)::

  >>> from sklearn.utils.estimator_checks import check_estimator
  >>> from sklearn.tree import DecisionTreeClassifier
  >>> check_estimator(DecisionTreeClassifier())  # passes
  [...]

The main motivation to make a class compatible to the scikit-learn estimator
interface might be that you want to use it together with model evaluation and
selection tools such as :class:`~model_selection.GridSearchCV` and
:class:`~pipeline.Pipeline`.

Before detailing the required interface below, we describe two ways to achieve
the correct interface more easily.

.. topic:: Project template:

    We provide a `project template
    <https://github.com/scikit-learn-contrib/project-template/>`_ which helps in the
    creation of Python packages containing scikit-learn compatible estimators. It
    provides:

    * an initial git repository with Python package directory structure
    * a template of a scikit-learn estimator
    * an initial test suite including use of :func:`~utils.parametrize_with_checks`
    * directory structures and scripts to compile documentation and example
      galleries
    * scripts to manage continuous integration (testing on Linux, MacOS, and Windows)
    * instructions from getting started to publishing on `PyPi <https://pypi.org/>`__

.. topic:: :class:`base.BaseEstimator` and mixins:

    We tend to use "duck typing" instead of checking for :func:`isinstance`, which means
    it's technically possible to implement estimator without inheriting from
    scikit-learn classes. However, if you don't inherit from the right mixins, either
    there will be a large amount of boilerplate code for you to implement and keep in
    sync with scikit-learn development, or your estimator might not function the same
    way as a scikit-learn estimator. Here we only document how to develop an estimator
    using our mixins. If you're interested in implementing your estimator without
    inheriting from scikit-learn mixins, you'd need to check our implementations.

    For example, below is a custom classifier, with more examples included in the
    scikit-learn-contrib `project template
    <https://github.com/scikit-learn-contrib/project-template/blob/master/skltemplate/_template.py>`__.

    It is particularly important to notice that mixins should be "on the left" while
    the ``BaseEstimator`` should be "on the right" in the inheritance list for proper
    MRO.

      >>> import numpy as np
      >>> from sklearn.base import BaseEstimator, ClassifierMixin
      >>> from sklearn.utils.validation import validate_data, check_is_fitted
      >>> from sklearn.utils.multiclass import unique_labels
      >>> from sklearn.metrics import euclidean_distances
      >>> class TemplateClassifier(ClassifierMixin, BaseEstimator):
      ...
      ...     def __init__(self, demo_param='demo'):
      ...         self.demo_param = demo_param
      ...
      ...     def fit(self, X, y):
      ...
      ...         # Check that X and y have correct shape, set n_features_in_, etc.
      ...         X, y = validate_data(self, X, y)
      ...         # Store the classes seen during fit
      ...         self.classes_ = unique_labels(y)
      ...
      ...         self.X_ = X
      ...         self.y_ = y
      ...         # Return the classifier
      ...         return self
      ...
      ...     def predict(self, X):
      ...
      ...         # Check if fit has been called
      ...         check_is_fitted(self)
      ...
      ...         # Input validation
      ...         X = validate_data(self, X, reset=False)
      ...
      ...         closest = np.argmin(euclidean_distances(X, self.X_), axis=1)
      ...         return self.y_[closest]

And you can check that the above estimator passes all common checks::

    >>> from sklearn.utils.estimator_checks import check_estimator
    >>> check_estimator(TemplateClassifier())  # passes            # doctest: +SKIP


get_params and set_params
-------------------------
All scikit-learn estimators have ``get_params`` and ``set_params`` functions.

The ``get_params`` function takes no arguments and returns a dict of the
``__init__`` parameters of the estimator, together with their values.

It takes one keyword argument, ``deep``, which receives a boolean value that determines
whether the method should return the parameters of sub-estimators (only relevant for
meta-estimators). The default value for ``deep`` is ``True``. For instance considering
the following estimator::

    >>> from sklearn.base import BaseEstimator
    >>> from sklearn.linear_model import LogisticRegression
    >>> class MyEstimator(BaseEstimator):
    ...     def __init__(self, subestimator=None, my_extra_param="random"):
    ...         self.subestimator = subestimator
    ...         self.my_extra_param = my_extra_param

The parameter `deep` controls control whether or not the parameters of the
`subestimator` should be reported. Thus when `deep=True`, the output will be::

    >>> my_estimator = MyEstimator(subestimator=LogisticRegression())
    >>> for param, value in my_estimator.get_params(deep=True).items():
    ...     print(f"{param} -> {value}")
    my_extra_param -> random
    subestimator__C -> 1.0
    subestimator__class_weight -> None
    subestimator__dual -> False
    subestimator__fit_intercept -> True
    subestimator__intercept_scaling -> 1
    subestimator__l1_ratio -> None
    subestimator__max_iter -> 100
    subestimator__multi_class -> deprecated
    subestimator__n_jobs -> None
    subestimator__penalty -> l2
    subestimator__random_state -> None
    subestimator__solver -> lbfgs
    subestimator__tol -> 0.0001
    subestimator__verbose -> 0
    subestimator__warm_start -> False
    subestimator -> LogisticRegression()

If the meta-estimator takes multiple sub-estimators, often, those sub-estimators have
names (as e.g. named steps in a :class:`~pipeline.Pipeline` object), in which case the
key should become `<name>__C`, `<name>__class_weight`, etc.

When ``deep=False``, the output will be::

    >>> for param, value in my_estimator.get_params(deep=False).items():
    ...     print(f"{param} -> {value}")
    my_extra_param -> random
    subestimator -> LogisticRegression()

On the other hand, ``set_params`` takes the parameters of ``__init__`` as keyword
arguments, unpacks them into a dict of the form ``'parameter': value`` and sets the
parameters of the estimator using this dict. It returns the estimator itself.

The :func:`~base.BaseEstimator.set_params` function is used to set parameters during
grid search for instance.

.. _cloning:

Cloning
-------
As already mentioned that when constructor arguments are mutable, they should be
copied before modifying them. This also applies to constructor arguments which are
estimators. That's why meta-estimators such as :class:`~model_selection.GridSearchCV`
create a copy of the given estimator before modifying it.

However, in scikit-learn, when we copy an estimator, we get an unfitted estimator
where only the constructor arguments are copied (with some exceptions, e.g. attributes
related to certain internal machinery such as metadata routing).

The function responsible for this behavior is :func:`~base.clone`.

Estimators can customize the behavior of :func:`base.clone` by overriding the
:func:`base.BaseEstimator.__sklearn_clone__` method. `__sklearn_clone__` must return an
instance of the estimator. `__sklearn_clone__` is useful when an estimator needs to hold
on to some state when :func:`base.clone` is called on the estimator. For example,
:class:`~sklearn.frozen.FrozenEstimator` makes use of this.

Estimator types
---------------
Among simple estimators (as opposed to meta-estimators), the most common types are
transformers, classifiers, regressors, and clustering algorithms.

**Transformers** inherit from :class:`~base.TransformerMixin`, and implement a `transform`
method. These are estimators which take the input, and transform it in some way. Note
that they should never change the number of input samples, and the output of `transform`
should correspond to its input samples in the same given order.

**Regressors** inherit from :class:`~base.RegressorMixin`, and implement a `predict` method.
They should accept numerical ``y`` in their `fit` method. Regressors use
:func:`~metrics.r2_score` by default in their :func:`~base.RegressorMixin.score` method.

**Classifiers** inherit from :class:`~base.ClassifierMixin`. If it applies, classifiers can
implement ``decision_function`` to return raw decision values, based on which
``predict`` can make its decision. If calculating probabilities is supported,
classifiers can also implement ``predict_proba`` and ``predict_log_proba``.

Classifiers should accept ``y`` (target) arguments to ``fit`` that are sequences (lists,
arrays) of either strings or integers. They should not assume that the class labels are
a contiguous range of integers; instead, they should store a list of classes in a
``classes_`` attribute or property. The order of class labels in this attribute should
match the order in which ``predict_proba``, ``predict_log_proba`` and
``decision_function`` return their values. The easiest way to achieve this is to put::

    self.classes_, y = np.unique(y, return_inverse=True)

in ``fit``.  This returns a new ``y`` that contains class indexes, rather than labels,
in the range [0, ``n_classes``).

A classifier's ``predict`` method should return arrays containing class labels from
``classes_``. In a classifier that implements ``decision_function``, this can be
achieved with::

    def predict(self, X):
        D = self.decision_function(X)
        return self.classes_[np.argmax(D, axis=1)]

The :mod:`~sklearn.utils.multiclass` module contains useful functions for working with
multiclass and multilabel problems.

**Clustering algorithms** inherit from :class:`~base.ClusterMixin`. Ideally, they should
accept a ``y`` parameter in their ``fit`` method, but it should be ignored. Clustering
algorithms should set a ``labels_`` attribute, storing the labels assigned to each
sample. If applicale, they can also implement a ``predict`` method, returning the
labels assigned to newly given samples.

If one needs to check the type of a given estimator, e.g. in a meta-estimator, one can
check if the given object implements a ``transform`` method for transformers, and
otherwise use helper functions such as :func:`~base.is_classifier` or
:func:`~base.is_regressor`.

.. _estimator_tags:

Estimator Tags
--------------
.. note::

    Scikit-learn introduced estimator tags in version 0.21 as a private API and mostly
    used in tests. However, these tags expanded over time and many third party
    developers also need to use them. Therefore in version 1.6 the API for the tags were
    revamped and exposed as public API.

The estimator tags are annotations of estimators that allow programmatic inspection of
their capabilities, such as sparse matrix support, supported output types and supported
methods. The estimator tags are an instance of :class:`~sklearn.utils.Tags` returned by
the method :meth:`~sklearn.base.BaseEstimator.__sklearn_tags__()`. These tags are used
in different places, such as :func:`~base.is_regressor` or the common checks run by
:func:`~sklearn.utils.estimator_checks.check_estimator` and
:func:`~sklearn.utils.estimator_checks.parametrize_with_checks`, where tags determine
which checks to run and what input data is appropriate. Tags can depend on estimator
parameters or even system architecture and can in general only be determined at runtime
and are therefore instance attributes rather than class attributes. See
:class:`~sklearn.utils.Tags` for more information about individual tags.

It is unlikely that the default values for each tag will suit the needs of your specific
estimator. You can change the default values by defining a `__sklearn_tags__()` method
which returns the new values for your estimator's tags. For example::

    class MyMultiOutputEstimator(BaseEstimator):

        def __sklearn_tags__(self):
            tags = super().__sklearn_tags__()
            tags.target_tags.single_output = False
            tags.non_deterministic = True
            return tags

You can create a new subclass of :class:`~sklearn.utils.Tags` if you wish to add new
tags to the existing set. Note that all attributes that you add in a child class need
to have a default value. It can be of the form::

    from dataclasses import dataclass, asdict

    @dataclass
    class MyTags(Tags):
        my_tag: bool = True

    class MyEstimator(BaseEstimator):
        def __sklearn_tags__(self):
            tags_orig = super().__sklearn_tags__()
            as_dict = {
                field.name: getattr(tags_orig, field.name)
                for field in fields(tags_orig)
            }
            tags = MyTags(**as_dict)
            tags.my_tag = True
            return tags


.. _developer_api_set_output:

Developer API for `set_output`
==============================

With
`SLEP018 <https://scikit-learn-enhancement-proposals.readthedocs.io/en/latest/slep018/proposal.html>`__,
scikit-learn introduces the `set_output` API for configuring transformers to
output pandas DataFrames. The `set_output` API is automatically defined if the
transformer defines :term:`get_feature_names_out` and subclasses
:class:`base.TransformerMixin`. :term:`get_feature_names_out` is used to get the
column names of pandas output.

:class:`base.OneToOneFeatureMixin` and
:class:`base.ClassNamePrefixFeaturesOutMixin` are helpful mixins for defining
:term:`get_feature_names_out`. :class:`base.OneToOneFeatureMixin` is useful when
the transformer has a one-to-one correspondence between input features and output
features, such as :class:`~preprocessing.StandardScaler`.
:class:`base.ClassNamePrefixFeaturesOutMixin` is useful when the transformer
needs to generate its own feature names out, such as :class:`~decomposition.PCA`.

You can opt-out of the `set_output` API by setting `auto_wrap_output_keys=None`
when defining a custom subclass::

    class MyTransformer(TransformerMixin, BaseEstimator, auto_wrap_output_keys=None):

        def fit(self, X, y=None):
            return self
        def transform(self, X, y=None):
            return X
        def get_feature_names_out(self, input_features=None):
            ...

The default value for `auto_wrap_output_keys` is `("transform",)`, which automatically
wraps `fit_transform` and `transform`. The `TransformerMixin` uses the
`__init_subclass__` mechanism to consume `auto_wrap_output_keys` and pass all other
keyword arguments to it's super class. Super classes' `__init_subclass__` should
**not** depend on `auto_wrap_output_keys`.

For transformers that return multiple arrays in `transform`, auto wrapping will
only wrap the first array and not alter the other arrays.

See :ref:`sphx_glr_auto_examples_miscellaneous_plot_set_output.py`
for an example on how to use the API.

.. _developer_api_check_is_fitted:

Developer API for `check_is_fitted`
===================================

By default :func:`~sklearn.utils.validation.check_is_fitted` checks if there
are any attributes in the instance with a trailing underscore, e.g. `coef_`.
An estimator can change the behavior by implementing a `__sklearn_is_fitted__`
method taking no input and returning a boolean. If this method exists,
:func:`~sklearn.utils.validation.check_is_fitted` simply returns its output.

See :ref:`sphx_glr_auto_examples_developing_estimators_sklearn_is_fitted.py`
for an example on how to use the API.

Developer API for HTML representation
=====================================

.. warning::

    The HTML representation API is experimental and the API is subject to change.

Estimators inheriting from :class:`~sklearn.base.BaseEstimator` display
a HTML representation of themselves in interactive programming
environments such as Jupyter notebooks. For instance, we can display this HTML
diagram::

    from sklearn.base import BaseEstimator

    BaseEstimator()

The raw HTML representation is obtained by invoking the function
:func:`~sklearn.utils.estimator_html_repr` on an estimator instance.

To customize the URL linking to an estimator's documentation (i.e. when clicking on the
"?" icon), override the `_doc_link_module` and `_doc_link_template` attributes. In
addition, you can provide a `_doc_link_url_param_generator` method. Set
`_doc_link_module` to the name of the (top level) module that contains your estimator.
If the value does not match the top level module name, the HTML representation will not
contain a link to the documentation. For scikit-learn estimators this is set to
`"sklearn"`.

The `_doc_link_template` is used to construct the final URL. By default, it can contain
two variables: `estimator_module` (the full name of the module containing the estimator)
and `estimator_name` (the class name of the estimator). If you need more variables you
should implement the `_doc_link_url_param_generator` method which should return a
dictionary of the variables and their values. This dictionary will be used to render the
`_doc_link_template`.

.. _coding-guidelines:

Coding guidelines
=================

The following are some guidelines on how new code should be written for
inclusion in scikit-learn, and which may be appropriate to adopt in external
projects. Of course, there are special cases and there will be exceptions to
these rules. However, following these rules when submitting new code makes
the review easier so new code can be integrated in less time.

Uniformly formatted code makes it easier to share code ownership. The
scikit-learn project tries to closely follow the official Python guidelines
detailed in `PEP8 <https://www.python.org/dev/peps/pep-0008>`_ that
detail how code should be formatted and indented. Please read it and
follow it.

In addition, we add the following guidelines:

* Use underscores to separate words in non class names: ``n_samples``
  rather than ``nsamples``.

* Avoid multiple statements on one line. Prefer a line return after
  a control flow statement (``if``/``for``).

* Use relative imports for references inside scikit-learn.

* Unit tests are an exception to the previous rule;
  they should use absolute imports, exactly as client code would.
  A corollary is that, if ``sklearn.foo`` exports a class or function
  that is implemented in ``sklearn.foo.bar.baz``,
  the test should import it from ``sklearn.foo``.

* **Please don't use** ``import *`` **in any case**. It is considered harmful
  by the `official Python recommendations
  <https://docs.python.org/3.1/howto/doanddont.html#at-module-level>`_.
  It makes the code harder to read as the origin of symbols is no
  longer explicitly referenced, but most important, it prevents
  using a static analysis tool like `pyflakes
  <https://divmod.readthedocs.io/en/latest/products/pyflakes.html>`_ to automatically
  find bugs in scikit-learn.

* Use the `numpy docstring standard
  <https://numpydoc.readthedocs.io/en/latest/format.html#docstring-standard>`_
  in all your docstrings.


A good example of code that we like can be found `here
<https://gist.github.com/nateGeorge/5455d2c57fb33c1ae04706f2dc4fee01>`_.

Input validation
----------------

.. currentmodule:: sklearn.utils

The module :mod:`sklearn.utils` contains various functions for doing input
validation and conversion. Sometimes, ``np.asarray`` suffices for validation;
do *not* use ``np.asanyarray`` or ``np.atleast_2d``, since those let NumPy's
``np.matrix`` through, which has a different API
(e.g., ``*`` means dot product on ``np.matrix``,
but Hadamard product on ``np.ndarray``).

In other cases, be sure to call :func:`check_array` on any array-like argument
passed to a scikit-learn API function. The exact parameters to use depends
mainly on whether and which ``scipy.sparse`` matrices must be accepted.

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

    from sklearn.utils import check_array, check_random_state

    def choose_random_sample(X, random_state=0):
        """Choose a random point from X.

        Parameters
        ----------
        X : array-like of shape (n_samples, n_features)
            An array representing the data.
        random_state : int or RandomState instance, default=0
            The seed of the pseudo random number generator that selects a
            random sample. Pass an int for reproducible output across multiple
            function calls.
            See :term:`Glossary <random_state>`.

        Returns
        -------
        x : ndarray of shape (n_features,)
            A random point selected from X.
        """
        X = check_array(X)
        random_state = check_random_state(random_state)
        i = random_state.randint(X.shape[0])
        return X[i]

If you use randomness in an estimator instead of a freestanding function,
some additional guidelines apply.

First off, the estimator should take a ``random_state`` argument to its
``__init__`` with a default value of ``None``.
It should store that argument's value, **unmodified**,
in an attribute ``random_state``.
``fit`` can call ``check_random_state`` on that attribute
to get an actual random number generator.
If, for some reason, randomness is needed after ``fit``,
the RNG should be stored in an attribute ``random_state_``.
The following example should make this clear::

    class GaussianNoise(BaseEstimator, TransformerMixin):
        """This estimator ignores its input and returns random Gaussian noise.

        It also does not adhere to all scikit-learn conventions,
        but showcases how to handle randomness.
        """

        def __init__(self, n_components=100, random_state=None):
            self.random_state = random_state
            self.n_components = n_components

        # the arguments are ignored anyway, so we make them optional
        def fit(self, X=None, y=None):
            self.random_state_ = check_random_state(self.random_state)

        def transform(self, X):
            n_samples = X.shape[0]
            return self.random_state_.randn(n_samples, self.n_components)

The reason for this setup is reproducibility:
when an estimator is ``fit`` twice to the same data,
it should produce an identical model both times,
hence the validation in ``fit``, not ``__init__``.

Numerical assertions in tests
-----------------------------

When asserting the quasi-equality of arrays of continuous values,
do use `sklearn.utils._testing.assert_allclose`.

The relative tolerance is automatically inferred from the provided arrays
dtypes (for float32 and float64 dtypes in particular) but you can override
via ``rtol``.

When comparing arrays of zero-elements, please do provide a non-zero value for
the absolute tolerance via ``atol``.

For more information, please refer to the docstring of
`sklearn.utils._testing.assert_allclose`.
