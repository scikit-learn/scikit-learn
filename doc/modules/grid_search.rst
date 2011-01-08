===================================================
Grid Search
===================================================

.. currentmodule:: scikits.learn.grid_search

`scikits.learn.grid_search` is a package to optimize
the parameters of a model (e.g. Support Vector Classifier)
using cross-validation.

The computation can be run in parallel using the multiprocessing package.

Main class is :class:`GridSearchCV`.

Examples
--------

See :ref:`example_grid_search_digits.py` for an example of
Grid Search computation on the digits dataset.

See :ref:`example_grid_search_text_feature_extraction.py` for an example
of Grid Search coupling parameters from a text documents feature extractor
(n-gram count vectorizer and TF-IDF transformer) with a classifier
(here a linear SVM trained with SGD with either elastic net or L2 penalty).

