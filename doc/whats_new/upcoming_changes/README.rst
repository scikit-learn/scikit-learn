:orphan:

Changelog
=========

This directory contains "news fragments" which are short files that contain a
small **ReST**-formatted text that will be added to the next what's new page.

Each file should be named like ``<PULL REQUEST>.<TYPE>.rst``, where
``<PULL REQUEST>`` is a pull request number, and ``<TYPE>`` is one of:

* ``major-feature``
* ``feature``
* ``efficiency``
* ``enhancement``
* ``fix``
* ``api``

See `this
<https://github.com/scikit-learn/scikit-learn/blob/main/doc/whats_new/changelog_legend.inc>`_
for more details about the meaning of each type.

If the change concerns a sub-package, the file should go in the sub-directory
relative to this sub-package, for example ``sklearn.linear_model`` or
``sklearn.tree``. In the vast majority of cases, your fragment should be
formatted as a bullet point.

So for example: ``28268.feature.rst`` would be added to the
``sklearn.ensemble`` folder and would have the content::

    - :class:`ensemble.ExtraTreesClassifier` and
      :class:`ensemble.ExtraTreesRegressor` now supports missing values in the
      data matrix `X`. Missing-values are handled by randomly moving all of the
      samples to the left, or right child node as the tree is traversed.
      :user:`Adam Li <adam2392>`

If you are unsure what pull request type to use, don't hesitate to ask in your
PR.

You can install ``towncrier`` and run ``towncrier --draft --version 1.6``
if you want to get a preview of how your change will look in the final release
notes.

There are a few other folders like ``changed-models`` ("Changes affecting many models" section), ``security``
("Security" section), etc ...

Type ``other`` is mostly meant to be used in the ``custom-top-level`` section.
