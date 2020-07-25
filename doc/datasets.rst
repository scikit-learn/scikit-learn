.. Places parent toc into the sidebar

:parenttoc: True

.. include:: includes/big_toc_css.rst

.. _datasets:

=========================
Dataset loading utilities
=========================

.. currentmodule:: sklearn.datasets

The ``sklearn.datasets`` package embeds some small toy datasets
as introduced in the :ref:`Getting Started <loading_example_dataset>` section.

This package also features helpers to fetch larger datasets commonly
used by the machine learning community to benchmark algorithms on data
that comes from the 'real world'.

To evaluate the impact of the scale of the dataset (``n_samples`` and
``n_features``) while controlling the statistical properties of the data
(typically the correlation and informativeness of the features), it is
also possible to generate synthetic data.

.. toctree::
    :maxdepth: 2

    datasets/general
    datasets/toy_dataset
    datasets/real_world
    datasets/sample_generators
    datasets/loading_other_datasets
