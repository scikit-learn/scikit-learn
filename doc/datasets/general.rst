.. Places parent toc into the sidebar

:parenttoc: True

General dataset API
===================

.. currentmodule:: sklearn.datasets

There are three main kinds of dataset interfaces that can be used to get
datasets depending on the desired type of dataset.

**The dataset loaders.** They can be used to load small standard datasets,
described in the :ref:`toy_datasets` section.

**The dataset fetchers.** They can be used to download and load larger datasets,
described in the :ref:`real_world_datasets` section.

Both loaders and fetchers functions return a :class:`~sklearn.utils.Bunch`
object holding at least two items:
an array of shape ``n_samples`` * ``n_features`` with
key ``data`` (except for 20newsgroups) and a numpy array of
length ``n_samples``, containing the target values, with key ``target``.

The Bunch object is a dictionary that exposes its keys are attributes.
For more information about Bunch object, see :class:`~sklearn.utils.Bunch`:

It's also possible for almost all of these function to constrain the output
to be a tuple containing only the data and the target, by setting the
``return_X_y`` parameter to ``True``.

The datasets also contain a full description in their ``DESCR`` attribute and
some contain ``feature_names`` and ``target_names``. See the dataset
descriptions below for details.

**The dataset generation functions.** They can be used to generate controlled
synthetic datasets, described in the :ref:`sample_generators` section.

These functions return a tuple ``(X, y)`` consisting of a ``n_samples`` *
``n_features`` numpy array ``X`` and an array of length ``n_samples``
containing the targets ``y``.

In addition, there are also miscellaneous tools to load datasets of other
formats or from other locations, described in the :ref:`loading_other_datasets`
section.
