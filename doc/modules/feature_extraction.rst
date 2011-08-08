
.. _feature_extraction:

==================
Feature extraction
==================

.. currentmodule:: scikits.learn.feature_extraction

The :mod:`scikits.learn.feature_extraction` module can be used to extract
features in a format supported by machine learning algorithms from datasets
consisting of formats such as text and image.


Text feature extraction
=======================

.. currentmodule:: scikits.learn.feature_extraction.text

XXX: a lot to do here


Image feature extraction
========================

.. currentmodule:: scikits.learn.feature_extraction.image

Patch extraction
----------------
The :func:`extract_patches_2d` and :class:`PatchExtractor` tools are useful for
extracting patches from two-dimensional images. For rebuilding an image from
all its patches, use :func:`reconstruct_from_patches_2d`.
