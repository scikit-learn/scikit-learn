..
    For doctests:

    >>> import numpy as np
    >>> import os
    >>> import tempfile
    >>> # Create a temporary folder for the data fetcher
    >>> custom_data_home = tempfile.mkdtemp()
    >>> os.makedirs(os.path.join(custom_data_home, 'mldata'))


.. _mldata:

Downloading datasets from the openml.org repository
===================================================

`openml.org <https://openml.org>`_ is a public repository for machine learning
data and experiments.

The ``sklearn.datasets`` package is able to directly download data
sets from the repository using the function
:func:`sklearn.datasets.fetch_openml`.

For example, to download a dataset of gene expressions in mice brains::

  >>> from sklearn.datasets import fetch_mldata
  >>> mice = fetch_mldata('miceprotein', data_home=custom_data_home)

The dataset contains a total of 70000 examples of handwritten digits
of size 28x28 pixels, labeled from 0 to 9::

  >>> mice.data.shape
  (70000, 784)
  >>> mice.target.shape
  (70000,)
  >>> np.unique(mice.target)
  array([ 0.,  1.,  2.,  3.,  4.,  5.,  6.,  7.,  8.,  9.])

You can get more information on the dataset by looking at the ``DESCR``
and ``details`` attributes::

  >>> print(mice.DESCR)
  something
  >>> mice.details

The ``DESCR`` contains a free-text description of the data, while ``details``
contains a dictionary of meta-data stored by openml, like the dataset id.
The id of the mice protein dataset is 4550, and you can use this (or the name)
to get more information on the dataset on the openml website: https://www.openml.org/d/4550.

Data sets in `mldata.org <http://mldata.org>`_ do not adhere to a strict
naming or formatting convention. :func:`sklearn.datasets.fetch_mldata` is
able to make sense of the most common cases, but allows to tailor the
defaults to individual datasets:

The id is also the best way to specify how to fetch a dataset from OpenML::

  >>> mice = fetch_mldata(4550, data_home=custom_data_home)
  >>> mice.details

Dataset Versions
----------------

A dataset is uniquely specified by its id, but not necessarily by its name.
Several different "versions" of a dataset with the same name can exist.
If a particular version of a dataset has been found to contain significant
issues, it might be inactivated. Using a name to specify a dataset will yield
the earliest version of a dataset that is still active. That means that
``fetch_mldata("miceprotein")`` can yield different results at differnt times
if earlier versions become inactive.
You can see that the dataset with id 4550 that we fetched above is the version 1
of the "miceprotein" dataset::

  >>> mice.details['version']
  1

In fact, this dataset only has one version. The iris dataset on the other hand
has multiple versions::

  >>> iris = fetch_mldata("iris", data_home=custom_data_home)
  >>> iris.details['version']
  >>> iris.details['id']

  >>> iris_61 = fetch_mldata(61, data_home=custom_data_home)
  >>> iris_61.details['version']
  >>> iris_61.details['id']

  >>> iris_969 = fetch_mldata(969, data_home=custom_data_home)
  >>> iris_969.details['version']
  >>> iris_969.details['id']

Specifying the dataset by the name "iris" yields the lowest version, version 1, with the id 61.
To make sure you always get this exact dataset, it is safest to specify it by the dataset id.
The other dataset, with id 969, is version 3 (version 2 has become inactive), and contains
a binarized version of the data::

  >>> np.bincount(iris_969.target)

..
    >>> import shutil
    >>> shutil.rmtree(custom_data_home)
