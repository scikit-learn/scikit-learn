#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
=========================================================
The Digit Dataset
=========================================================
This dataset is made up of 1797 8x8 images. Each image,
like the one shown below, is of a hand-written digit.
In order to utilize an 8x8 figure like this, we'd have to
first transform it into a feature vector with length 64.

See `here
<http://archive.ics.uci.edu/ml/datasets/Pen-Based+Recognition+of+Handwritten+Digits>`_
for more information about this dataset.
"""
print(__doc__)


# Code source: Gael Varoqueux
# Modified for Documentation merge by Jaques Grobler
# License: BSD 3 clause

from sklearn import datasets

import pylab as pl

#Load the digits dataset
digits = datasets.load_digits()

#Display the last digit
pl.figure(1, figsize=(3, 3))
pl.imshow(digits.images[-1], cmap=pl.cm.gray_r, interpolation='nearest')
pl.show()

from sklearn.datasets import make_translated_images

#Translates the last digits, show the "1px up" translation
X, y = make_translated_images([digits.data[-1]],
                              [digits.target[-1]],
                              images_shape=(8, 8), n_samples=4)
pl.imshow(X[-1].reshape((8, 8)), cmap=pl.cm.gray_r, interpolation='nearest')
pl.show()
