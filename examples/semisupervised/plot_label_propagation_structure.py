"""
===================================================
Plot Label Propagation learning a complex structure
===================================================

Example of LabelPropagation learning a complex internal structure
to demonstrate "manifold learning". The outer circle should be
labeled "red" and the inner circle "blue". Because both label groups
lie inside their own distinct shape, we can see that the labels
propagate correctly around the circle.
"""
print __doc__

from sklearn import label_propagation
import numpy as np
import pylab as pl

from sklearn import svm, datasets

label_spread = label_propagation.LabelSpreading(gamma=20, alpha=1.0)

# generate ring with inner box
n_samples_per_circle = 100
outer_circ_xs = np.cos(np.linspace(0, 2 * np.pi, n_samples_per_circle))
outer_circ_ys = np.sin(np.linspace(0, 2 * np.pi, n_samples_per_circle))
inner_circ_xs = np.cos(np.linspace(0, 2 * np.pi, n_samples_per_circle)) * 0.45
inner_circ_ys = np.sin(np.linspace(0, 2 * np.pi, n_samples_per_circle)) * 0.45

all_xs = np.append(outer_circ_xs, inner_circ_xs)
all_ys = np.append(outer_circ_ys, inner_circ_ys)
data = np.vstack((np.append(outer_circ_xs, inner_circ_xs),\
       np.append(outer_circ_ys, inner_circ_ys))).T
labels = ['outer'] +\
         ['unlabeled' for x in range(0, n_samples_per_circle - 1)] +\
         ['inner'] +\
         ['unlabeled' for x in range(0, n_samples_per_circle - 1)]

label_spread.fit(data, labels, unlabeled_identifier='unlabeled')

output_labels = label_spread.transduction_

pl.subplot(2, 1, 1)
plot_outer_labeled, = pl.plot(outer_circ_xs[0], outer_circ_ys[0], 'rs-')
plot_unlabeled, = pl.plot(np.append(outer_circ_xs[1:], inner_circ_xs[1:]), \
        np.append(outer_circ_ys[1:], inner_circ_ys[1:]), 'g.')
plot_inner_labeled, = pl.plot(inner_circ_xs[0], inner_circ_ys[0], 'bs-')
pl.legend((plot_outer_labeled, plot_inner_labeled, plot_unlabeled), \
        ('Outer Labeled', 'Inner Labeled', 'Unlabeled'), 'upper left', \
        numpoints=1, shadow=False)
pl.title("Raw data (2 classes=red and blue)")

pl.subplot(2, 1, 2)
output_label_array = np.asarray(output_labels)
outer_numbers = np.where(output_label_array == 'outer')
inner_numbers = np.where(output_label_array == 'inner')
plot_outer, = pl.plot(all_xs[outer_numbers], all_ys[outer_numbers], 'rs-')
plot_inner, = pl.plot(all_xs[inner_numbers], all_ys[inner_numbers], 'bs-')
pl.legend((plot_outer, plot_inner), ('Outer Learned', 'Inner Learned'), \
        'upper left', numpoints=1, shadow=False)
pl.title("Labels learned with Label Spreading")

pl.show()
