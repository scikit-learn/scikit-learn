.. _density_estimation:

==================
Density Estimation
==================
.. sectionauthor:: Jake Vanderplas <vanderplas@astro.washington.edu>

Density estimation walks the line between unsupervised learning, feature
engineering, and data modeling.  Some of the most popular and useful
density estimation techniques are mixture models such as
Gaussian Mixtures (:class:`~sklearn.mixture.GaussianMixture`), and
neighbor-based approaches such as the kernel density estimate
(:class:`~sklearn.neighbors.KernelDensity`).
Gaussian Mixtures are discussed more fully in the context of
:ref:`clustering <clustering>`, because the technique is also useful as
an unsupervised clustering scheme.

Density estimation is a very simple concept, and most people are already
familiar with one common density estimation technique: the histogram.

Density Estimation: Histograms
==============================
A histogram is a simple visualization of data where bins are defined, and the
number of data points within each bin is tallied.  An example of a histogram
can be seen in the upper-left panel of the following figure:

.. |hist_to_kde| image:: ../auto_examples/neighbors/images/sphx_glr_plot_kde_1d_001.png
   :target: ../auto_examples/neighbors/plot_kde_1d.html
   :scale: 80

.. centered:: |hist_to_kde|

A major problem with histograms, however, is that the choice of binning can
have a disproportionate effect on the resulting visualization.  Consider the
upper-right panel of the above figure.  It shows a histogram over the same
data, with the bins shifted right.  The results of the two visualizations look
entirely different, and might lead to different interpretations of the data.

Intuitively, one can also think of a histogram as a stack of blocks, one block
per point.  By stacking the blocks in the appropriate grid space, we recover
the histogram.  But what if, instead of stacking the blocks on a regular grid,
we center each block on the point it represents, and sum the total height at
each location?  This idea leads to the lower-left visualization.  It is perhaps
not as clean as a histogram, but the fact that the data drive the block
locations mean that it is a much better representation of the underlying
data.

This visualization is an example of a *kernel density estimation*, in this case
with a top-hat kernel (i.e. a square block at each point).  We can recover a
smoother distribution by using a smoother kernel.  The bottom-right plot shows
a Gaussian kernel density estimate, in which each point contributes a Gaussian
curve to the total.  The result is a smooth density estimate which is derived
from the data, and functions as a powerful non-parametric model of the
distribution of points.

.. _kernel_density:

Kernel Density Estimation
=========================
Kernel density estimation in scikit-learn is implemented in the
:class:`~sklearn.neighbors.KernelDensity` estimator, which uses the
Ball Tree or KD Tree for efficient queries (see :ref:`neighbors` for
a discussion of these).  Though the above example
uses a 1D data set for simplicity, kernel density estimation can be
performed in any number of dimensions, though in practice the curse of
dimensionality causes its performance to degrade in high dimensions.

In the following figure, 100 points are drawn from a bimodal distribution,
and the kernel density estimates are shown for three choices of kernels:

.. |kde_1d_distribution| image:: ../auto_examples/neighbors/images/sphx_glr_plot_kde_1d_003.png
   :target: ../auto_examples/neighbors/plot_kde_1d.html
   :scale: 80

.. centered:: |kde_1d_distribution|

It's clear how the kernel shape affects the smoothness of the resulting
distribution.  The scikit-learn kernel density estimator can be used as
follows:

   >>> from sklearn.neighbors import KernelDensity
   >>> import numpy as np
   >>> X = np.array([[-1, -1], [-2, -1], [-3, -2], [1, 1], [2, 1], [3, 2]])
   >>> kde = KernelDensity(kernel='gaussian', bandwidth=0.2).fit(X)
   >>> kde.score_samples(X)
   array([-0.41075698, -0.41075698, -0.41076071, -0.41075698, -0.41075698,
          -0.41076071])

Here we have used ``kernel='gaussian'``, as seen above.
Mathematically, a kernel is a positive function :math:`K(x;h)`
which is controlled by the bandwidth parameter :math:`h`.
Given this kernel form, the density estimate at a point :math:`y` within
a group of points :math:`x_i; i=1\cdots N` is given by:

.. math::
    \rho_K(y) = \sum_{i=1}^{N} K(y - x_i; h)

The bandwidth here acts as a smoothing parameter, controlling the tradeoff
between bias and variance in the result.  A large bandwidth leads to a very
smooth (i.e. high-bias) density distribution.  A small bandwidth leads
to an unsmooth (i.e. high-variance) density distribution.

:class:`~sklearn.neighbors.KernelDensity` implements several common kernel
forms, which are shown in the following figure:

.. |kde_kernels| image:: ../auto_examples/neighbors/images/sphx_glr_plot_kde_1d_002.png
   :target: ../auto_examples/neighbors/plot_kde_1d.html
   :scale: 80

.. centered:: |kde_kernels|

The form of these kernels is as follows:

* Gaussian kernel (``kernel = 'gaussian'``)

  :math:`K(x; h) \propto \exp(- \frac{x^2}{2h^2} )`

* Tophat kernel (``kernel = 'tophat'``)

  :math:`K(x; h) \propto 1` if :math:`x < h`

* Epanechnikov kernel (``kernel = 'epanechnikov'``)

  :math:`K(x; h) \propto 1 - \frac{x^2}{h^2}`

* Exponential kernel (``kernel = 'exponential'``)

  :math:`K(x; h) \propto \exp(-x/h)`

* Linear kernel (``kernel = 'linear'``)

  :math:`K(x; h) \propto 1 - x/h` if :math:`x < h`

* Cosine kernel (``kernel = 'cosine'``)

  :math:`K(x; h) \propto \cos(\frac{\pi x}{2h})` if :math:`x < h`

The kernel density estimator can be used with any of the valid distance
metrics (see :class:`~sklearn.neighbors.DistanceMetric` for a list of available metrics), though
the results are properly normalized only for the Euclidean metric.  One
particularly useful metric is the
`Haversine distance <https://en.wikipedia.org/wiki/Haversine_formula>`_
which measures the angular distance between points on a sphere.  Here
is an example of using a kernel density estimate for a visualization
of geospatial data, in this case the distribution of observations of two
different species on the South American continent:

.. |species_kde| image:: ../auto_examples/neighbors/images/sphx_glr_plot_species_kde_001.png
   :target: ../auto_examples/neighbors/plot_species_kde.html
   :scale: 80

.. centered:: |species_kde|

One other useful application of kernel density estimation is to learn a
non-parametric generative model of a dataset in order to efficiently
draw new samples from this generative model.
Here is an example of using this process to
create a new set of hand-written digits, using a Gaussian kernel learned
on a PCA projection of the data:

.. |digits_kde| image:: ../auto_examples/neighbors/images/sphx_glr_plot_digits_kde_sampling_001.png
   :target: ../auto_examples/neighbors/plot_digits_kde_sampling.html
   :scale: 80

.. centered:: |digits_kde|

The "new" data consists of linear combinations of the input data, with weights
probabilistically drawn given the KDE model.

.. topic:: Examples:

  * :ref:`sphx_glr_auto_examples_neighbors_plot_kde_1d.py`: computation of simple kernel
    density estimates in one dimension.

  * :ref:`sphx_glr_auto_examples_neighbors_plot_digits_kde_sampling.py`: an example of using
    Kernel Density estimation to learn a generative model of the hand-written
    digits data, and drawing new samples from this model.

  * :ref:`sphx_glr_auto_examples_neighbors_plot_species_kde.py`: an example of Kernel Density
    estimation using the Haversine distance metric to visualize geospatial data
