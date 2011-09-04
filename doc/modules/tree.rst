
.. _tree:

==============
Decision Trees
==============

.. currentmodule:: sklearn.tree

**Decision Trees** are a supervised learning
method used for :ref:`classification <tree_classification>`,
:ref:`regression <tree_regression>`.

The advantages of Decision Trees are:

    - Simple to understand and interpret.  Trees can be visualised if needed.
    
    - Requires little data preparation. Other techniques often require data normalisation, dummy variables need to be created and blank values to be removed.
    
    - Able to handle both numerical and categorical data. Other techniques are usually specialised in analysing datasets that have only one type of variable.
    
    - Uses a white box model. If a given situation is observable in a model the explanation for the condition is easily explained by boolean logic. An example of a black box model is an artificial neural network since the explanation for the results is difficult to understand.
    
    - Possible to validate a model using statistical tests. That makes it possible to account for the reliability of the model.
    
    - Robust. Performs well even if its assumptions are somewhat violated by the true model from which the data were generated.

The disadvantages of Decision Trees include:

    - The problem of learning an optimal decision tree is known to be NP-complete under several aspects of optimality and even for simple concepts. Consequently, practical decision-tree learning algorithms are based on heuristic algorithms such as the greedy algorithm where locally optimal decisions are made at each node. Such algorithms cannot guarantee to return the globally optimal decision tree.
    
    - Decision-tree learners can create over-complex trees that do not generalise the data well. This is called overfitting. Mechanisms such as pruning are necessary to avoid this problem.
    
    - There are concepts that are hard to learn because decision trees do not express them easily, such as XOR, parity or multiplexer problems. 

    - Decision tree learners create biased trees if some classes dominate. It is therefore recommended to balance the dataset prior to fitting with the Decision Tree

.. _tree_classification:

Classification
==============

:class:`DecisionTreeClassifier` is capable of performing multi-class classification on a dataset.


.. figure:: ../auto_examples/tree/images/plot_iris_1.png
   :target: ../auto_examples/tree/plot_iris.html
   :align: center


    >>> from sklearn import tree
    >>> X = [[0, 0], [1, 1]]
    >>> Y = [0, 1]
    >>> clf = tree.DecisionTreeClassifier()
    >>> clf = clf.fit(X, Y)


After being fitted, the model can then be used to predict new values::

    >>> clf.predict([[2., 2.]])
    array([1])

Unbalanced problems
--------------------



.. topic:: Examples:

 * :ref:`example_tree_plot_iris.py`,
 * :ref:`example_tree_plot_separating_hyperplane.py`,


.. _tree_regression:

Regression
==========

Decision Trees can be applied to solve regression problems. 

Uses Mean Squared Error to find the best fit.

As with classification classes, the fit method will take as
argument vectors X, y, only that in this case y is expected to have
floating point values instead of integer values::

    >>> from sklearn import tree
    >>> X = [[0, 0], [2, 2]]
    >>> y = [0.5, 2.5]
    >>> clf = tree.DecisionTreeRegressor()
    >>> clf = clf.fit(X, y)
    >>> clf.predict([[1, 1]])
    array([ 2.5])


.. topic:: Examples:

 * :ref:`example_tree_plot_tree_regression.py`


Tips on Practical Use
=====================

.. _tree_mathematical_formulation:

Mathematical formulation
========================

.. topic:: References:
