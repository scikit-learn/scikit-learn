=========================
Dataset loading utilities
=========================

.. currentmodule:: scikits.learn.datasets

The ``scikits.learn.datasets`` package embeds some small toy datasets
as introduced in the "Getting Started" section.

To evaluate the impact of the scale of the dataset (``n_features`` and
``n_samples``) while controlling the statistical properties of the data
(typically the correlation and informativeness of the features), it is
also possible to generate synthetic data data

This package also features helpers to fetch larger datasets commonly
used by the machine learning community to benchmark algorithm on data
that comes from the 'real world'.


Dataset generators
==================

TODO



The Labeled Faces in the Wild face recognition dataset
======================================================

This dataset is a collection of JPEG pictures of famous people collected
over the internet, all details are available on the official website:

    http://vis-www.cs.umass.edu/lfw/

Each picture is centered on a single face. The typical task is called
Face Verification: given a pair of two pictures, a binary classifier
must predict whether the two images are from the same person.

An alternative task, Face Recognition or Face Identification is:
given the picture of the face of an unknown person, identify the name
of the person by refering to a gallery of previously seen pictures of
identified persons.

Both Face Verification and Face Recognition are tasks that are typically
performed on the output of a model trained to perform Face Detection. The
most popular model for Face Detection is called Viola-Johns and is
implemented in the OpenCV library. The LFW faces were extracted by this
face detector from various online websites.

``scikit-learn`` provides two loaders that will automatically download,
cache, parse the metadata files, decode the jpeg and convert the
interesting slices into memmaped numpy arrays. This dataset size if more
than 200 MB. The first load typically takes more than a couple of minutes
to fully decode the relevant part of the JPEG files into numpy arrays. If
the dataset has  been loaded once, the following times the loading times
less than 200ms by using a memmaped version memoized on the disk in the
``~/scikit_learn_data/lfw_home/`` folder using ``joblib``.

The first loader is used for the Face Identification task: a multi-class
classification task (hence supervised learning)::

  >>> from scikits.learn.datasets import load_lfw_people
  >>> lfw_people = load_lfw_people(min_faces_per_person=100)

  >>> for name in lfw_people.class_names:
  ...     print name
  ...
  Colin Powell
  Donald Rumsfeld
  George W Bush
  Gerhard Schroeder
  Tony Blair

The default slice is a rectangular shape around the face, removing
most of the background::

  >>> lfw_people.data.dtype
  dtype('float32')

  >>> lfw_people.data.shape
  (1140, 62, 47)

Each of the ``1140`` faces is assigned to a single person id in the ``target``
array::

  >>> lfw_people.target.shape
  (1140,)

  >>> lfw_people.target[:10]
  memmap([2, 3, 1, 4, 1, 0, 2, 0, 2, 1])

The second loader is typically used for the face verification task: each sample
is a pair of two picture belonging or not to the same person::

  >>> from scikits.learn.datasets import load_lfw_pairs
  >>> lfw_pairs_train = load_lfw_pairs(subset='train')

  >>> list(lfw_pairs_train.class_names)
  ['Different persons', 'Same person']

  >>> lfw_pairs_train.data.shape
  (2200, 2, 62, 47)

  >>> lfw_pairs_train.target.shape
  (2200,)

Both for the ``load_lfw_people`` and ``load_lfw_pairs`` function it is
possible to get an additional dimension with the RGB color channels by
passing ``color=True``::

  >>> lfw_pairs_train = load_lfw_pairs(subset='train', color=True)
  >>> lfw_pairs_train.data.shape
  (2200, 2, 62, 47, 3)

The ``load_lfw_pairs`` datasets is subdived in 3 subsets: the development
``train`` set, the development ``test`` set and an evaluation ``10_folds``
set meant to be compute performance metrics using a 10-folds cross
validation scheme.


Examples
--------

:ref:`example_applications_plot_face_recognition.py`


The 20 newsgroups text dataset
==============================

TODO




