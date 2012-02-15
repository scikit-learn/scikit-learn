# Code source: Gael Varoqueux
# Modified for Documentation merge by Jaques Grobler
# License: BSD

try:
    from sklearn import datasets
except ImportError:
    from scikits.learn import datasets

import pylab as pl

#Load the digits dataset
digits = datasets.load_digits()

#Display the first digit
pl.imshow(digits.images[0], cmap=pl.cm.gray_r)
pl.show()

