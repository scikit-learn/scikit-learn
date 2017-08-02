.. _freeze:

Frozen estimators and transfer learning
=======================================

.. currentmodule:: sklearn

It can be useful to pre-fit an estimator before including it in a Pipeline,
FeatureUnion or other meta-estimators.  Example applications include:

* transfer learning: incorporating a transformer trained on a large unlabelled
  dataset in a prediction pipeline where the data to be modelled is much smaller
* feature selection on the basis of an already fitted predictive model

To enable this, your estimator can be wrapped in :class:`freeze.FreezeWrap`.
For example::

    Without transfer learning

    >>> from sklearn.datasets import load_...
    >>> from sklearn.model_selection import cross_val_score
    >>> cross_val_score(make_pipeline(TfidfVectorizer(), LogisticRegression()),
    ...                 X, y)

    With transfer learning:
    >>> from sklearn.freeze import FreezeWrap
    >>> tfidf = TfidfVectorizer().fit(large_X)
    >>> cross_val_score(make_pipeline(FreezeWrap(tfidf), LogisticRegression()),
    ...                 X, y)

In particular, calling ``FrezeWrap(tfidf).fit(X, y)`` now does nothing,
while calling ``FrezeWrap(tfidf).fit_transform(X, y)`` just returns the result of
``tfidf.transform(X)``.

.. note::
    When an estimator is frozen, calling :func:`clone` on it will return
    itself.::

        >>> from base import clone
        >>> frozen = FreezeWrap(tfidf)
        >>> clone(frozen) is frozen
        True

    This allows the model to be left untouched in cross-validation and
    meta-estimators which clear the estimator with ``clone``.

.. warning:: Leakage:
    Please take care to not introduce data leakage by this method: do not
    incorporate your test set into the training of some frozen component,
    unless it would be realistic to do so in the target application.
