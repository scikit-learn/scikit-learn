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

For example, to download a dataset of gene expressions in mice brains:

  >>> from sklearn.datasets import fetch_mldata
  >>> mnist = fetch_mldata('miceprotein', data_home=custom_data_home)

The MNIST database contains a total of 70000 examples of handwritten digits
of size 28x28 pixels, labeled from 0 to 9::

  >>> mnist.data.shape
  (70000, 784)
  >>> mnist.target.shape
  (70000,)
  >>> np.unique(mnist.target)
  array([ 0.,  1.,  2.,  3.,  4.,  5.,  6.,  7.,  8.,  9.])

After the first download, the dataset is cached locally in the path
specified by the ``data_home`` keyword argument, which defaults to
``~/scikit_learn_data/``::

  >>> os.listdir(os.path.join(custom_data_home, 'mldata'))
  ['mnist-original.mat']

Data sets in `mldata.org <http://mldata.org>`_ do not adhere to a strict
naming or formatting convention. :func:`sklearn.datasets.fetch_mldata` is
able to make sense of the most common cases, but allows to tailor the
defaults to individual datasets:

* For datasets with multiple columns, :func:`sklearn.datasets.fetch_mldata`
  tries to identify the target and data columns and rename them to ``target``
  and ``data``. This is done by looking for arrays named ``label`` and
  ``data`` in the dataset, and failing that by choosing the first array to be
  ``target`` and the second to be ``data``. This behavior can be changed with
  the ``target_name`` and ``data_name`` keywords, setting them to a specific
  name or index number (the name and order of the columns in the datasets
  can be found at its `mldata.org <http://mldata.org>`_ under the tab "Data"::

    >>> iris2 = fetch_mldata('datasets-UCI iris', target_name=1, data_name=0,
    ...                      data_home=custom_data_home)
    >>> iris3 = fetch_mldata('datasets-UCI iris', target_name='class',
    ...                      data_name='double0', data_home=custom_data_home)


..
    >>> import shutil
    >>> shutil.rmtree(custom_data_home)
