# ruff: noqa: CPY001
"""
=======================================
Release Highlights for scikit-learn 1.8
=======================================

.. currentmodule:: sklearn

We are pleased to announce the release of scikit-learn 1.8! Many bug fixes
and improvements were added, as well as some key new features. Below we
detail the highlights of this release. **For an exhaustive list of
all the changes**, please refer to the :ref:`release notes <release_notes_1_8>`.

To install the latest version (with pip)::

    pip install --upgrade scikit-learn

or with conda::

    conda install -c conda-forge scikit-learn

"""

# %%
# Array API stuff
# ---------------
# TODO copy and paste from 1.7 highlights needs tweaking
# Several functions have been updated to support array API compatible inputs since
# version 1.7, especially TODO from the :mod:`sklearn.metrics` module.
#
# Please refer to the :ref:`array API support<array_api>` page for instructions to use
# scikit-learn with array API compatible libraries such as PyTorch or CuPy.
#
# TODO are there more "important ones"? could be in the example see point below
# TODO Only show highlighted code without executing it since we don't have a
# GPU in the doc build? We could also show snippet PyTorch CPU with
# commented out device='cuda' if you want to run on GPU you only have to
# uncomment it. Alternative idea link to Colab notebook?

# %%
# Free-threaded CPython 3.14 support
# ----------------------------------
#
# scikit-learn has support for free-threaded CPython, in particular
# free-threaded wheels are available for all of our supported platforms on Python
# 3.14.
#
# Free-threaded (also known as nogil) CPython is a version of CPython that aims at
# enabling efficient multi-threaded use cases by removing the Global Interpreter
# Lock (GIL).
#
# If you want to try out free-threaded Python, the recommendation is to use
# Python 3.14, that has fixed a number of issues compared to Python 3.13. Feel
# free to try free-threaded on your use case and report any issues!
#
# For more details about free-threaded CPython see `py-free-threading doc <https://py-free-threading.github.io>`_,
# in particular `how to install a free-threaded CPython <https://py-free-threading.github.io/installing_cpython/>`_
# and `Ecosystem compatibility tracking <https://py-free-threading.github.io/tracking/>`_.

# %%
# Temperature scaling in `CalibratedClassifierCV`
# -----------------------------------------------
# TODO

# %%
# Linear models improvements
# --------------------------
# TODO

# %%
# HTML representation of estimators
# ---------------------------------
# TODO

# %%
# DecisionTreeRegressor with MAE
#

# %%
# ClassicalMDS
# ------------
