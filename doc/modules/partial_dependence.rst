
.. _partial_dependence:

==================
Partial dependence
==================

.. currentmodule:: sklearn.inspection

Partial dependence plots (PDP) and individual conditional
expectation (ICE) plots show the dependence between the target
response [1]_ and a set of 'target' features, marginalizing over the values
of all other features (the 'complement' features). While PDPs show the
average effect of the 'target' features, ICE plots visualize the dependence
of the prediction on the 'target' features for each instance separately,
resulting in one curve per instance. Intuitively, we can
interpret the partial dependence as the expected target response as a
function of the 'target' features.

Due to the limits of human perception the size of the target feature set
must be small (usually, one for ICE and one or two for PDP) thus the target
features are usually chosen among the most important features.

The figure below shows four one-way and one two-way partial dependence plots
for the California housing dataset, with a :class:`GradientBoostingRegressor
<sklearn.ensemble.GradientBoostingRegressor>`:

.. figure:: ../auto_examples/inspection/images/sphx_glr_plot_partial_dependence_002.png
   :target: ../auto_examples/inspection/plot_partial_dependence.html
   :align: center
   :scale: 70

One-way PDPs tell us about the interaction between the target response and
the target feature (e.g. linear, non-linear). The upper left plot in the
above figure shows the effect of the median income in a district on the
median house price; we can clearly see a linear relationship among them. Note
that PDPs assume that the target features are independent from the complement
features, and this assumption is often violated in practice.

PDPs with two target features show the interactions among the two features.
For example, the two-variable PDP in the above figure shows the dependence
of median house price on joint values of house age and average occupants per
household. We can clearly see an interaction between the two features: for
an average occupancy greater than two, the house price is nearly independent of
the house age, whereas for values less than 2 there is a strong dependence
on age.

The figure below shows the individual conditional expectation plots for the
same dataset, with a :class:`MLPRegressor
<sklearn.neural_network.MLPRegressor>`:

.. figure:: ../auto_examples/inspection/images/sphx_glr_plot_individual_conditional_expectation_002.png
   :target: ../auto_examples/inspection/plot_individual_conditional_expectation.html
   :align: center
   :scale: 70

While the PDPs show the average effect of the target features, ICE plots
effect of the target variable for each instance separately. For example,
we could observe a linear relationship among median income and the house
price in PDP. However, the corresponding ICE plot show that there are some
exceptions, where the house price remain constant with the median income.
Similar to PDPs, the target features are assumed to be independent from the
complement features.

The :mod:`sklearn.inspection` module provides a convenience functions,
:func:`plot_partial_dependence` and
:func:`plot_individual_conditional_expectation` to create partial
dependence plots and individual conditional expectation plots respectively.

In the below example we show how to create a grid of
partial dependence plots: two one-way PDPs for the features ``0`` and ``1``
and a two-way PDP between the two features::

    >>> from sklearn.datasets import make_hastie_10_2
    >>> from sklearn.ensemble import GradientBoostingClassifier
    >>> from sklearn.inspection import plot_partial_dependence

    >>> X, y = make_hastie_10_2(random_state=0)
    >>> clf = GradientBoostingClassifier(n_estimators=100, learning_rate=1.0,
    ...     max_depth=1, random_state=0).fit(X, y)
    >>> features = [0, 1, (0, 1)]
    >>> plot_partial_dependence(clf, X, features) #doctest: +SKIP

Individual conditional expectation plots for the features ``0`` and ``1`` can
be generated with::

    >>> from sklearn.inspection import plot_individual_conditional_expectation
    >>> features = [0, 1]
    >>> plot_individual_conditional_expectation(clf, X, features) #doctest: +SKIP

You can access the newly created figure and Axes objects using ``plt.gcf()``
and ``plt.gca()``.

For multi-class classification, you need to set the class label for which
the PDP/ICEs should be created via the ``target`` argument::

    >>> from sklearn.datasets import load_iris
    >>> iris = load_iris()
    >>> mc_clf = GradientBoostingClassifier(n_estimators=10,
    ...     max_depth=1).fit(iris.data, iris.target)
    >>> features = [3, 2, (3, 2)]
    >>> plot_partial_dependence(mc_clf, X, features, target=0) #doctest: +SKIP
    >>> features = [3, 2]
    >>> plot_individual_conditional_expectation(mc_clf, X, features, target=0) #doctest: +SKIP

The same parameter ``target`` is used to specify the target in multi-output
regression settings.

If you need the raw values of the partial dependence/individual conditional
expectation function rather than the plots, you can use functions
:func:`sklearn.inspection.partial_dependence` and
:func:`individual_conditional_expectation`::

    >>> from sklearn.inspection import partial_dependence

    >>> pdp, axes = partial_dependence(clf, X, [0])
    >>> pdp
    array([[ 2.466...,  2.466..., ...
    >>> axes
    [array([-1.624..., -1.592..., ...

The values at which the partial dependence/individual conditional
expectation should be evaluated are directly generated from ``X``.
For 2-way partial dependence, a 2D-grid of values is generated.
The ``values`` field returned by :func:`sklearn.inspection.partial_dependence`
or :func:`individual_conditional_expectation` gives the actual values
used in the grid for each target feature. They also correspond to the axis
of the plots.

Mathematical Definition
^^^^^^^^^^^^^^^^^^^^^^^

Let :math:`X_S` be the set of target features (i.e. the `features` parameter)
and let :math:`X_C` be its complement.

The partial dependence of the response :math:`f` at a point :math:`x_S` is
defined as:

.. math::

    pd_{X_S}(x_S) &\overset{def}{=} \mathbb{E}_{X_C}\left[ f(x_S, X_C) \right]\\
                  &= \int f(x_S, x_C) p(x_C) dx_C,

where :math:`f(x_S, x_C)` is the response function (:term:`predict`,
:term:`predict_proba` or :term:`decision_function`) for a given sample whose
values are defined by :math:`x_S` for the features in :math:`X_S`, and by
:math:`x_C` for the features in :math:`X_C`. Note that :math:`x_S` and
:math:`x_C` may be tuples.

Computing this integral for various values of :math:`x_S` produces a partial
dependence plot as above.

In ICE, for each instance in :math:`\{(x_S^{(i)}, x_C^{(i)})\}_{i=1}^N` the
curve :math:`f_S^{(i)}` is calculated against :math:`x_S^{(i)}`, while
:math:`x_C^{(i)}` remains fixed.

Computation methods
^^^^^^^^^^^^^^^^^^^

There are two main methods to approximate the integral above, namely the
'brute' and 'recursion' methods. The `method` parameter controls which method
to use.

The 'brute' method is a generic method that works with any estimator. It
approximates the above integral by computing an average over the data `X`:

.. math::

    pd_{X_S}(x_S) \approx \frac{1}{n_\text{samples}} \sum_{i=1}^n f(x_S, x_C^{(i)}),

where :math:`x_C^{(i)}` is the value of the i-th sample for the features in
:math:`X_C`. For each value of :math:`x_S`, this method requires a full pass
over the dataset `X` which is computationally intensive. Note that only
'brute' method is available for individual conditional expectation
calculation, since the averaging is already implicit in the 'recursion' method.

The 'recursion' method is faster than the 'brute' method, but it is only
supported by some tree-based estimators. It is computed as follows. For a
given point :math:`x_S`, a weighted tree traversal is performed: if a split
node involves a 'target' feature, the corresponding left or right branch is
followed; otherwise both branches are followed, each branch being weighted
by the fraction of training samples that entered that branch. Finally, the
partial dependence is given by a weighted average of all the visited leaves
values.

With the 'brute' method, the parameter `X` is used both for generating the
grid of values :math:`x_S` and the complement feature values :math:`x_C`.
However with the 'recursion' method, `X` is only used for the grid values:
implicitly, the :math:`x_C` values are those of the training data.

By default, the 'recursion' method is used on tree-based estimators that
support it, and 'brute' is used for the rest.

.. _pdp_method_differences:

.. note::

    While both methods should be close in general, they might differ in some
    specific settings. The 'brute' method assumes the existence of the
    data points :math:`(x_S, x_C^{(i)})`. When the features are correlated,
    such artificial samples may have a very low probability mass. The 'brute'
    and 'recursion' methods will likely disagree regarding the value of the
    partial dependence, because they will treat these unlikely
    samples differently. Remember, however, that the primary assumption for
    interpreting PDPs is that the features should be independent.

.. rubric:: Footnotes

.. [1] For classification, the target response may be the probability of a
   class (the positive class for binary classification), or the decision
   function.

.. topic:: Examples:

 * :ref:`sphx_glr_auto_examples_inspection_plot_partial_dependence.py`
 * :ref:`sphx_glr_auto_examples_inspection_plot_individual_conditional_expectation.py`

.. topic:: References

    T. Hastie, R. Tibshirani and J. Friedman, `The Elements of
    Statistical Learning <https://web.stanford.edu/~hastie/ElemStatLearn//>`_,
    Second Edition, Section 10.13.2, Springer, 2009.

    C. Molnar, `Interpretable Machine Learning
    <https://christophm.github.io/interpretable-ml-book/>`_, Section 5.1, 2019.
