..
    For doctests:

    >>> import numpy as np
    >>> import os
    >>> import tempfile
    >>> # Create a temporary folder for the data fetcher
    >>> custom_data_home = tempfile.mkdtemp()
    >>> os.makedirs(os.path.join(custom_data_home, 'openml'))


.. _openml:

Downloading datasets from the openml.org repository
===================================================

`openml.org <https://openml.org>`_ is a public repository for machine learning
data and experiments, that allows everybody to upload open datasets.

The ``sklearn.datasets`` package is able to directly download data
sets from the repository using the function
:func:`sklearn.datasets.fetch_openml`.

For example, to download a dataset of gene expressions in mice brains::

  >>> from sklearn.datasets import fetch_openml
  >>> mice = fetch_openml('miceprotein', data_home=custom_data_home)

The dataset contains a total of 70000 examples of handwritten digits
of size 28x28 pixels, labeled from 0 to 9::

  >>> mice.data.shape
  (1080, 81)
  >>> mice.target.shape
  (1080,)
  >>> np.unique(mice.target) # doctest: +NORMALIZE_WHITESPACE +ELLIPSIS
  array([b"'c-CS-m'", b"'c-CS-s'", b"'c-SC-m'", b"'c-SC-s'", b"'t-CS-m'",
  b"'t-CS-s'", b"'t-SC-m'", b"'t-SC-s'"], dtype='|S8')

You can get more information on the dataset by looking at the ``DESCR``
and ``details`` attributes::

  >>> print(mice.DESCR) # doctest: +NORMALIZE_WHITESPACE +ELLIPSIS
  **Author**: Clara Higuera, Katheleen J. Gardiner, Krzysztof J. Cios  
  **Source**: [UCI](https://archive.ics.uci.edu/ml/datasets/Mice+Protein+Expression) - 2015   
  **Please cite**: Higuera C, Gardiner KJ, Cios KJ (2015) Self-Organizing
  Feature Maps Identify Proteins Critical to Learning in a Mouse Model of Down
  Syndrome. PLoS ONE 10(6): e0129126...

  >>> mice.details # doctest: +NORMALIZE_WHITESPACE +ELLIPSIS
  {'id': '4550', 'name': 'MiceProtein', 'version': '1', 'format': 'ARFF',
  'creator': ...,
  'upload_date': '2016-02-17T14:32:49', 'licence': 'Public', 'url':
  'https://www.openml.org/data/v1/download/1804243/MiceProtein.ARFF', 'file_id':
  '1804243', 'default_target_attribute': 'class', 'citation': 'Higuera C,
  Gardiner KJ, Cios KJ (2015) Self-Organizing Feature Maps Identify Proteins
  Critical to Learning in a Mouse Model of Down Syndrome. PLoS ONE 10(6):
  e0129126. [Web Link] journal.pone.0129126', 'tag': ['OpenML100', 'study_14',
  'study_34'], 'visibility': 'public', 'status': 'active', 'md5_checksum':
  '3c479a6885bfa0438971388283a1ce32'}


The ``DESCR`` contains a free-text description of the data, while ``details``
contains a dictionary of meta-data stored by openml, like the dataset id.
The id of the mice protein dataset is 4550, and you can use this (or the name)
to get more information on the dataset on the openml website: https://www.openml.org/d/4550.

The id is also the best way to specify how to fetch a dataset from OpenML::

  >>> mice = fetch_openml(4550, data_home=custom_data_home)
  >>> mice.details # doctest: +NORMALIZE_WHITESPACE +ELLIPSIS
  {'id': '4550', 'name': 'MiceProtein', 'version': '1', 'format': 'ARFF',
  'creator': ...,
  'upload_date': '2016-02-17T14:32:49', 'licence': 'Public', 'url':
  'https://www.openml.org/data/v1/download/1804243/MiceProtein.ARFF', 'file_id':
  '1804243', 'default_target_attribute': 'class', 'citation': 'Higuera C,
  Gardiner KJ, Cios KJ (2015) Self-Organizing Feature Maps Identify Proteins
  Critical to Learning in a Mouse Model of Down Syndrome. PLoS ONE 10(6):
  e0129126. [Web Link] journal.pone.0129126', 'tag': ['OpenML100', 'study_14',
  'study_34'], 'visibility': 'public', 'status': 'active', 'md5_checksum':
  '3c479a6885bfa0438971388283a1ce32'}

Dataset Versions
----------------

A dataset is uniquely specified by its id, but not necessarily by its name.
Several different "versions" of a dataset with the same name can exist.
If a particular version of a dataset has been found to contain significant
issues, it might be inactivated. Using a name to specify a dataset will yield
the earliest version of a dataset that is still active. That means that
``fetch_openml("miceprotein")`` can yield different results at differnt times
if earlier versions become inactive.
You can see that the dataset with id 4550 that we fetched above is the version 1
of the "miceprotein" dataset::

  >>> mice.details['version']
  '1'

In fact, this dataset only has one version. The iris dataset on the other hand
has multiple versions::

  >>> iris = fetch_openml("iris", data_home=custom_data_home)
  >>> iris.details['version']
  '1'
  >>> iris.details['id']
  '61'

  >>> iris_61 = fetch_openml(61, data_home=custom_data_home)
  >>> iris_61.details['version']
  '1'
  >>> iris_61.details['id']
  '61'

  >>> iris_969 = fetch_openml(969, data_home=custom_data_home)
  >>> iris_969.details['version']
  '3'
  >>> iris_969.details['id']
  '969'

'Specifying the dataset by the name "iris" yields the lowest version, version 1, with the id 61.
To make sure you always get this exact dataset, it is safest to specify it by the dataset id.
The other dataset, with id 969, is version 3 (version 2 has become inactive), and contains
a binarized version of the data::

  >>> np.unique(iris_969.target)
  array([b'N', b'P'],
        dtype='|S1')

..
    >>> import shutil
    >>> shutil.rmtree(custom_data_home)
