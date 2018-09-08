..
    For doctests:

    >>> import numpy as np
    >>> import os


.. _openml:

Downloading datasets from the openml.org repository
===================================================

`openml.org <https://openml.org>`_ is a public repository for machine learning
data and experiments, that allows everybody to upload open datasets.

The ``sklearn.datasets`` package is able to download datasets
from the repository using the function
:func:`sklearn.datasets.fetch_openml`.

For example, to download a dataset of gene expressions in mice brains::

  >>> from sklearn.datasets import fetch_openml
  >>> mice = fetch_openml(name='miceprotein', version=4)

To fully specify a dataset, you need to provide a name and a version, though
the version is optional, see :ref:`openml_versions` below.
The dataset contains a total of 1080 examples belonging to 8 different
classes::

  >>> mice.data.shape
  (1080, 77)
  >>> mice.target.shape
  (1080,)
  >>> np.unique(mice.target) # doctest: +NORMALIZE_WHITESPACE
  array(['c-CS-m', 'c-CS-s', 'c-SC-m', 'c-SC-s', 't-CS-m', 't-CS-s', 't-SC-m', 't-SC-s'], dtype=object)

You can get more information on the dataset by looking at the ``DESCR``
and ``details`` attributes::

  >>> print(mice.DESCR) # doctest: +NORMALIZE_WHITESPACE +ELLIPSIS +SKIP
  **Author**: Clara Higuera, Katheleen J. Gardiner, Krzysztof J. Cios
  **Source**: [UCI](https://archive.ics.uci.edu/ml/datasets/Mice+Protein+Expression) - 2015
  **Please cite**: Higuera C, Gardiner KJ, Cios KJ (2015) Self-Organizing
  Feature Maps Identify Proteins Critical to Learning in a Mouse Model of Down
  Syndrome. PLoS ONE 10(6): e0129126...

  >>> mice.details # doctest: +NORMALIZE_WHITESPACE +ELLIPSIS +SKIP
  {'id': '40966', 'name': 'MiceProtein', 'version': '4', 'format': 'ARFF',
  'upload_date': '2017-11-08T16:00:15', 'licence': 'Public',
  'url': 'https://www.openml.org/data/v1/download/17928620/MiceProtein.arff',
  'file_id': '17928620', 'default_target_attribute': 'class',
  'row_id_attribute': 'MouseID',
  'ignore_attribute': ['Genotype', 'Treatment', 'Behavior'],
  'tag': ['OpenML-CC18', 'study_135', 'study_98', 'study_99'],
  'visibility': 'public', 'status': 'active',
  'md5_checksum': '3c479a6885bfa0438971388283a1ce32'}


The ``DESCR`` contains a free-text description of the data, while ``details``
contains a dictionary of meta-data stored by openml, like the dataset id.
For more details, see the `OpenML documentation
<https://docs.openml.org/#data>`_ The ``data_id`` of the mice protein dataset
is 40966, and you can use this (or the name) to get more information on the
dataset on the openml website::

  >>> mice.url
  'https://www.openml.org/d/40966'

The ``data_id`` also uniquely identifies a dataset from OpenML::

  >>> mice = fetch_openml(data_id=40966)
  >>> mice.details # doctest: +NORMALIZE_WHITESPACE +ELLIPSIS +SKIP
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

.. _openml_versions:

Dataset Versions
----------------

A dataset is uniquely specified by its ``data_id``, but not necessarily by its
name. Several different "versions" of a dataset with the same name can exist
which can contain entirely different datasets.
If a particular version of a dataset has been found to contain significant
issues, it might be deactivated. Using a name to specify a dataset will yield
the earliest version of a dataset that is still active. That means that
``fetch_openml(name="miceprotein")`` can yield different results at different
times if earlier versions become inactive.
You can see that the dataset with ``data_id`` 40966 that we fetched above is
the version 1 of the "miceprotein" dataset::

  >>> mice.details['version']  #doctest: +SKIP
  '1'

In fact, this dataset only has one version. The iris dataset on the other hand
has multiple versions::

  >>> iris = fetch_openml(name="iris")
  >>> iris.details['version']  #doctest: +SKIP
  '1'
  >>> iris.details['id']  #doctest: +SKIP
  '61'

  >>> iris_61 = fetch_openml(data_id=61)
  >>> iris_61.details['version']
  '1'
  >>> iris_61.details['id']
  '61'

  >>> iris_969 = fetch_openml(data_id=969)
  >>> iris_969.details['version']
  '3'
  >>> iris_969.details['id']
  '969'

Specifying the dataset by the name "iris" yields the lowest version, version 1,
with the ``data_id`` 61. To make sure you always get this exact dataset, it is
safest to specify it by the dataset ``data_id``. The other dataset, with
``data_id`` 969, is version 3 (version 2 has become inactive), and contains a
binarized version of the data::

  >>> np.unique(iris_969.target)
  array(['N', 'P'], dtype=object)

You can also specify both the name and the version, which also uniquely
identifies the dataset::

  >>> iris_version_3 = fetch_openml(name="iris", version=3)
  >>> iris_version_3.details['version']
  '3'
  >>> iris_version_3.details['id']
  '969'


.. topic:: References:

 * Vanschoren, van Rijn, Bischl and Torgo
   `"OpenML: networked science in machine learning"
   <https://arxiv.org/pdf/1407.7722.pdf>`_,
   ACM SIGKDD Explorations Newsletter, 15(2), 49-60, 2014.
