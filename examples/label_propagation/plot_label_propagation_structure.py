"""
===================================================
Plot Label Propagation learning a complex structure
===================================================

Example of LabelPropagation learning a complex internal structure
"""
print __doc__

import numpy as np
import pylab as pl
from scikits.learn import svm, datasets

X = iris.data
Y = iris.target

lp = label_propagation.LabelPropagation()
lspread = label_propagation.LabelSpreading()

# title for the plots
titles = [ ]

#TODO ALL OF THIS
