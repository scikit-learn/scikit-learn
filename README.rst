.. -*- mode: rst -*-

|Azure|_ |CirrusCI|_ |Codecov|_ |CircleCI|_ |Nightly wheels|_ |Black|_ |PythonVersion|_ |PyPi|_ |DOI|_ |Benchmark|_

.. |Azure| image:: https://dev.azure.com/scikit-learn/scikit-learn/_apis/build/status/scikit-learn.scikit-learn?branchName=main
.. _Azure: https://dev.azure.com/scikit-learn/scikit-learn/_build/latest?definitionId=1&branchName=main

.. |CircleCI| image:: https://circleci.com/gh/scikit-learn/scikit-learn/tree/main.svg?style=shield&circle-token=:circle-token
.. _CircleCI: https://circleci.com/gh/scikit-learn/scikit-learn

.. |CirrusCI| image:: https://img.shields.io/cirrus/github/scikit-learn/scikit-learn/main?label=Cirrus%20CI
.. _CirrusCI: https://cirrus-ci.com/github/scikit-learn/scikit-learn/main

.. |Codecov| image:: https://codecov.io/gh/scikit-learn/scikit-learn/branch/main/graph/badge.svg?token=Pk8G9gg3y9
.. _Codecov: https://codecov.io/gh/scikit-learn/scikit-learn

.. |Nightly wheels| image:: https://github.com/scikit-learn/scikit-learn/workflows/Wheel%20builder/badge.svg?event=schedule
.. _`Nightly wheels`: https://github.com/scikit-learn/scikit-learn/actions?query=workflow%3A%22Wheel+builder%22+event%3Aschedule

.. |PythonVersion| image:: https://img.shields.io/badge/python-3.8%20%7C%203.9%20%7C%203.10-blue
.. _PythonVersion: https://pypi.org/project/scikit-learn/

.. |PyPi| image:: https://img.shields.io/pypi/v/scikit-learn
.. _PyPi: https://pypi.org/project/scikit-learn

.. |Black| image:: https://img.shields.io/badge/code%20style-black-000000.svg
.. _Black: https://github.com/psf/black

.. |DOI| image:: https://zenodo.org/badge/21369/scikit-learn/scikit-learn.svg
.. _DOI: https://zenodo.org/badge/latestdoi/21369/scikit-learn/scikit-learn

.. |Benchmark| image:: https://img.shields.io/badge/Benchmarked%20by-asv-blue
.. _`Benchmark`: https://scikit-learn.org/scikit-learn-benchmarks/

.. |PythonMinVersion| replace:: 3.8
.. |NumPyMinVersion| replace:: 1.17.3
.. |SciPyMinVersion| replace:: 1.5.0
.. |JoblibMinVersion| replace:: 1.1.1
.. |ThreadpoolctlMinVersion| replace:: 2.0.0
.. |MatplotlibMinVersion| replace:: 3.1.3
.. |Scikit-ImageMinVersion| replace:: 0.16.2
.. |PandasMinVersion| replace:: 1.0.5
.. |SeabornMinVersion| replace:: 0.9.0
.. |PytestMinVersion| replace:: 5.3.1
.. |PlotlyMinVersion| replace:: 5.10.0

=================
Scikit-learn-tree
=================

``scikit-learn-tree`` is a maintained fork of scikit-learn, which advances the tree submodule, while staying in-line
with changes from upstream scikit-learn. It is an exact stand-in for ``sklearn`` in package imports, but is
released under the name ``scikit-learn-tree`` to avoid confusion.

It is currently maintained by a team of volunteers.

The upstream package **scikit-learn** is a Python module for machine learning built on top of
SciPy and is distributed under the 3-Clause BSD license. Refer to their website for all documentation
needs: https://scikit-learn.org.

Why a fork?
-----------
Currently, the scikit-learn tree submodule is difficult to extend. Requests to modularize
and improve the extensibility of the code is currently unsupported, or may take a long time.
The desire for advanced tree models that also leverage the robustness of scikit-learn is desirable.

However, "hard-forking" via copy/pasting the explicit Python/Cython code into another tree package
altogether is undesirable because it results in a tree codebase that is inherently different
and not compatible with ``scikit-learn``. For example, `quantile-forests <https://github.com/zillow/quantile-forest>`_,
and `EconML <https://github.com/py-why/EconML>`_ do this, and their current tree submodules
cannot take advantage of improvements made in upstream ``scikit-learn``.

An example of seamless integration would be `scikit-survival <https://github.com/sebp/scikit-survival>`_, which
only needs to implement a subclass of the Cython ``Criterion`` oject in their code to enable survival trees.

Maintaining a "soft-fork" of ``scikit-learn`` in the form of a repository fork allows us to develop
a separate package that serves as a stand-in for ``sklearn`` in any package, extends the tree submodule
and can also be synced with upstream changes in ``scikit-learn``. This enables this fork to always
take advantage of improvements made in ``scikit-learn`` main upstream, while providing a customizable
tree API.

Installation
------------

Dependencies
~~~~~~~~~~~~

scikit-learn requires:

- Python (>= |PythonMinVersion|)
- NumPy (>= |NumPyMinVersion|)
- SciPy (>= |SciPyMinVersion|)
- joblib (>= |JoblibMinVersion|)
- threadpoolctl (>= |ThreadpoolctlMinVersion|)

============================
Installing scikit-learn-tree
============================

Scikit-learn-tree is a maintained fork of scikit-learn, which extends the
tree submodule in a few ways documented in `fork_changelog`_. 

We release versions of scikit-learn-tree in an analagous fashion to
scikit-learn main. Due to maintenance resources, we only release on PyPi
and recommend therefore installing with ``pip``.

There are different ways to install scikit-learn-tree:

  * Install the latest official release `install_fork_release`_. This
    is the best approach for most users. It will provide a stable version
    and pre-built packages are available for most platforms.
    
  * Building the package from source `install_source`_. This is best for users who want the
    latest-and-greatest features and aren't afraid of running
    brand-new code. This is also needed for users who wish to contribute to the
    project.

.. _install_fork_release:

Installing the latest release
-----------------------------
We release wheels for common distributions and this is thus installable via pip.

    pip install scikit-learn-tree

This will install ``scikit-learn-tree`` under the namespace of ``sklearn``, which then
can be used as a stand-in for any package that relies on the public API of ``sklearn``.

For example, any usage of ``scikit-learn`` is preserved with ``scikit-learn-tree``

  >>> # the sklearn installed is that of scikit-learn-tree and is equivalent to scikit-learn
  >>> from sklearn.ensemble import RandomForestClassifier
  >>> clf = RandomForestClassifier(random_state=0)
  >>> X = [[ 1,  2,  3],  # 2 samples, 3 features
  ...      [11, 12, 13]]
  >>> y = [0, 1]  # classes of each sample
  >>> clf.fit(X, y)
  RandomForestClassifier(random_state=0)

.. _install_source:

Building from source
--------------------
If you are a developer and are interested in helping maintain, or add some new
features to the fork, the building from source instructions are exactly the same
as that of scikit-learn main, so please refer to `scikit-learn documentation <https://scikit-learn.org/stable/developers/advanced_installation.html#install-bleeding-edge>`_
for instructions on building from source.

===========

Development
-----------

We welcome new contributors of all experience levels, specifically to maintain the fork.
Any contributions that make sure our fork is "better in-line" with scikit-learn upstream,
or improves the tree submodule in anyway will be appreciated.

The scikit-learn community goals are to be helpful, welcoming, and effective. The
`Development Guide <https://scikit-learn.org/stable/developers/index.html>`_
has detailed information about contributing code, documentation, tests, and
more. We've included some basic information in this README.

=========================

.. _fork_changelog:

Major Changes of the Fork
-------------------------

The purpose of this page is to illustrate some of the main features that
``scikit-learn-tree`` provides compared to ``scikit-learn``. It assumes a
an understanding of core package ``scikit-learn`` and also decision trees
models. Please refer to our installation instructions `install_fork_release`_ for installing ``scikit-learn-tree``.

Scikit-learn-tree though operates as a stand-in for upstream ``scikit-learn``.
It is used in packages exactly the same way and will support all features
in the corresponding version of ``scikit-learn``. For example, if you
are interested in features of ``scikit-learn`` in v1.2.2 for ``NearestNeighbors`` algorithm,
then if ``scikit-learn-tree`` has a version release of v1.2.2, then it will have
all those features. 

The breaking API changes will be with respect to anything in the ``tree`` submodule,
and related Forest ensemble models. See below for a detailed list of breaking changes.

See: https://scikit-learn.org/ for documentation on scikit-learn main.

Our Philosophy
--------------
Our design philosophy with this fork of ``scikit-learn`` is to maintain as few changes
as possible, such that incorporating upstream changes into the fork requires minimal effort.

Candidate changes and PRs accepted into the fork are those that:

- improve compatability with upstream ``scikit-learn`` main
- enable improved extensibility of tree models

Decision tree generalizations
-----------------------------

``Scikit-learn`` provides an axis-aligned `sklearn.tree.DecisionTreeClassifier <https://scikit-learn.org/stable/modules/generated/sklearn.tree.DecisionTreeClassifier.html>`_
decision tree model (classifier and regressor), which has a few fundamental limitations
that prevent 3rd parties from utilizing the existing class, without forking a large
amount of copy/pasted Python and Cython code. We highlight those limitations here
and then describe how we generalize that limitation.

Cython Internal Private API:

Note, the Cython API for scikit-learn is still not a publicly supported API, so it may
change without warning.

- leaf and split nodes: These nodes are treated the same way and there is no internal
  API for setting them differently. Quantile trees and causal trees inherently generalize
  how leaf nodes are set.
- Criterion class: The criterion class currently assumes a supervised learning interface.
  - Our fix: We implement a ``BaseCriterion`` object that provides an abstract API for unsupervised criterion.
- Splitter class: The splitter clas currently assumes a supervised learning interface and
  does not provide a way of generalizing the way split candidates are proposed.
  - Our fix: We implement a ``BaseSplitter`` object that provides an abstract API for unsupervised splitters and also implement an API to allow generalizations of the ``SplitRecord`` struct and ``Splitter.node_split`` function. For example, this enables oblique splits to be considered.
- Tree class: The tree class currently assumes a supervised learning interface and does not
  provide a way of generalizing the type of tree.
  - Our fix: We implementa ``BaseTree`` object that provides an abstract API for general tree models and also implement an API that allows generalization of the type of tree. For example, oblique trees are trivially implementable as an extension now.
- stopping conditions for splitter: Currently, the ``Splitter.node_split`` function has various
  stopping conditions for the splitter based on hyperparameters. It is plausible that these conditions
  may be extended. For example, in causal trees, one may want the splitter to also account for
  a minimal degree of heterogeneity (i.e. variance) in its children nodes. 

Python API:

- ``sklearn.tree.BaseDecisionTree`` assumes the underlying tree model is supervised: The ``y``
  parameter is required to be passed in, which is not necessary for general tree-based models.
  For example, an unsupervised tree may pass in ``y=None``.
  - Our fix: We fix this API, so the ``BaseDecisionTree`` is subclassable by unsupervised tree models that do not require ``y`` to be defined.
- ``sklearn.tree.BaseDecisionTree`` does not provide a way to generalize the ``Criterion``, ``Splitter``
  and ``Tree`` Cython classes used: The current codebase requires users to define custom
  criterion and/or splitters outside the instantiation of the ``BaseDecisionTree``. This prevents
  users from generalizing the ``Criterion`` and ``Splitter`` and creating a neat Python API wrapper.
  Moreover, the ``Tree`` class is not customizable.
  - Our fix: We internally implement a private function to actually build the entire tree, ``BaseDecisionTree._build_tree``, which can be overridden in subclasses that customize the criterion, splitter, or tree, or any combination of them.
- ``sklearn.ensemble.BaseForest`` and its subclass algorithms are slow when ``n_samples`` is very high. Binning
  features into a histogram, which is the basis of "LightGBM" and "HistGradientBoostingClassifier" is a computational
  trick that can both significantly increase runtime efficiency, but also help prevent overfitting in trees, since
  the sorting in "BestSplitter" is done on bins rather than the continuous feature values. This would enable
  random forests and their variants to scale to millions of samples.
  - Our fix: We added a ``max_bins=None`` keyword argument to the ``BaseForest`` class, and all its subclasses. The default behavior is no binning. The current implementation is not necessarily efficient. There are several improvements to be made. See below.

Overall, the existing tree models, such as `sklearn.tree.DecisionTreeClassifier <https://scikit-learn.org/stable/modules/generated/sklearn.tree.DecisionTreeClassifier.html>`_
and `sklearn.ensemble.RandomForestClassifier <https://scikit-learn.org/stable/modules/generated/sklearn.ensemble.RandomForestClassifier.html#sklearn.ensemble.RandomForestClassifier>`_ all work exactly the same as they
would in ``scikit-learn`` main, but these extensions enable 3rd-party packages to extend
the Cython/Python API easily.

Roadmap
-------
There are several improvements that can be made in this fork. Primarily, the binning feature
promises to make Random Forests and their variants ultra-fast. However, the binning needs
to be implemented in a similar fashion to ``HistGradientBoostingClassifier``, which passes
in the binning thresholds throughout the tree construction step, such that the split nodes
store the actual numerical value of the bin rather than the "bin index". This requires
modifying the tree Cython code to take in a ``binning_thresholds`` parameter that is part
of the ``_BinMapper`` fitted class. This also allows us not to do any binning during prediction/apply
time because the tree already stores the "numerical" threshold value we would want to apply
to any incoming ``X`` that is not binned.

Besides that modification, the tree and splitter need to be able to handle not just ``np.float32``
data (the type for X normally in Random Forests), but also ``uint8`` data (the type for X when it
is binned in to e.g. 255 bins). This would not only save RAM since ``uint8`` storage of millions
of samples would result in many GB saved, but also improved runtime.

So in summary, the Cython code of the tree submodule needs to take in an extra parameter for
the binning thresholds if binning occurs and also be able to handle ``X`` being of dtype ``uint8``.
Afterwards, Random Forests will have fully leveraged the binning feature.

Something to keep in mind is that upstream scikit-learn is actively working on incorporating
missing-value handling and categorical handling into Random Forests.

Next steps
----------

We have briefly covered how the tree submodule has changed with respect to ``scikit-learn``.
This enables packages to leverage these changes in developing more complex tree models
that may, or may not eventually be PRed into ``scikit-learn``. For example,

- `scikit-tree <https://docs.neurodata.io/scikit-tree/dev/index.html>`_ is a scikit-learn
  compatible package for more complex and advanced tree models.

If you are developing tree models, we encourage you to take a look at that package, or
if you have suggestions to make the tree submodule of our fork, ``scikit-learn-tree``
more 