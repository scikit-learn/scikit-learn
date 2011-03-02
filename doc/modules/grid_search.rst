===================================================
Grid Search
===================================================

.. currentmodule:: scikits.learn.grid_search


Grid Search is used to optimize the parameters of a model
(e.g. Support Vector Classifier, Lasso, etc.) using cross-validation.

Main class is :class:`GridSearchCV`.

Examples
--------

See :ref:`example_grid_search_digits.py` for an example of
Grid Search computation on the digits dataset.

See :ref:`example_grid_search_text_feature_extraction.py` for an example
of Grid Search coupling parameters from a text documents feature extractor
(n-gram count vectorizer and TF-IDF transformer) with a classifier
(here a linear SVM trained with SGD with either elastic net or L2 penalty).

Notes
-----
Computations can be run in parallel if your OS supports it, by using
the keyword n_jobs=-1, see function signature for more details.
