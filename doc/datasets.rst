.. _datasets:

=========================
Dataset loading utilities
=========================

.. currentmodule:: sklearn.datasets

The ``sklearn.datasets`` package embeds some small toy datasets and provides helpers
to fetch larger datasets commonly used by the machine learning community to benchmark
algorithms on data that comes from the 'real world'.

To evaluate the impact of the scale of the dataset (``n_samples`` and
``n_features``) while controlling the statistical properties of the data
(typically the correlation and informativeness of the features), it is
also possible to generate synthetic data.

**General dataset API.** There are three main kinds of dataset interfaces that
can be used to get datasets depending on the desired type of dataset.

**The dataset loaders.** They can be used to load small standard datasets,
described in the :ref:`toy_datasets` section.

**The dataset fetchers.** They can be used to download and load larger datasets,
described in the :ref:`real_world_datasets` section.

Both loaders and fetchers functions return a :class:`~sklearn.utils.Bunch`
object holding at least two items:
an array of shape ``n_samples`` * ``n_features`` with
key ``data`` (except for 20newsgroups) and a numpy array of
length ``n_samples``, containing the target values, with key ``target``.

The Bunch object is a dictionary that exposes its keys as attributes.
For more information about Bunch object, see :class:`~sklearn.utils.Bunch`.

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


.. toctree::
    :maxdepth: 2

    datasets/toy_dataset
    datasets/real_world
    datasets/sample_generators
    datasets/loading_other_datasets
