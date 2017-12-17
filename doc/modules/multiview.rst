
.. currentmodule:: sklearn.multiview

.. _multiview:

==================
Multiview learning
==================


Multiview learning provides multiview methods to work with multiview data
(datasets with several data matrices from the same samples). It contains methods
for multiview dimensionality reduction and methods for multiview clustering.


.. figure:: ../auto_examples/multiview/images/sphx_glr_plot_mvmds_001.png
   :target: ../auto_examples/multiview/plot_mvmds.html
   :align: center
   :scale: 60

The multiview learning implementations available in scikit-learn are summarized below.


Introduction
============

Given a multiview dataset with ``v`` input data matrices, multiview dimensionality
reduction methods produce a single, low-dimensional projection of the input data 
samples, trying to maintain as much of the original information as possible.

Multiview learning can be thought of as an effort to adapt some of the manifold algorithms
principles to algorithms whose main inputs are multiview data. These data can be understood
as different views of the same data, so it is attempted to use these views to and all their
information for dimensionality reduction and spectral clustering.

Methods developed here are adaptions of single view algorithms to multiview data.
Also, these modules are translation from the multiview package, firstly written in
R.

.. _mvmds:

MVMDS
=====

MVMDS (Multiview Multidimensional Scaling) ia one of the approaches to dimensionality reduction that offers the class 
:class:`mvmds` to perform multiview dimensionality reduction in a similar way than
the multidimensional scaling method does (similar to ``cmdscale`` in R language).
In general, it is a technique used for analyzing similarity or dissimilarity data.

.. figure:: ../auto_examples/multiview/images/sphx_glr_plot_mvmds_001.png
   :target: ../auto_examples/multiview/plot_mvmds.html
   :align: center
   :scale: 60

Complexity
----------
MVMDS computes distance matrices from plain data matrices is so and preprocesses these distance
matrices (centering and double square). Lastly, it extracts commom principal components
of the processed data. The overall complexity is about
:math:`O[k n^2 v]`.

* :math:`n`: number of samples of each view
* :math:`k`: components to extract

.. _mvtsne:

MvtSNE
======

Another dimensionality reduction function in this package is the class :class:`mvtsne`, that
extends the ``tsne`` algorithm (available in ``manifold`` module) to work with multiview data.
It based on the conversion of affinities of data to probabilities. The affinities
in the original space are represented by Gaussian joint probabilities and the affinities
in the embedded space are represented by Student’s t-distributions.

.. figure:: ../auto_examples/multiview/images/sphx_glr_plot_mvtsne_001.png
   :target: ../auto_examples/multiview/plot_mvtsne.html
   :align: center
   :scale: 60

Complexity
----------
The MvtSNE algorithm comprises two steps:

1. **Opinion pooling finding**: MvtSNE computes optimal pooling of probabilities for
each view from affinities of data points. It runs an optimizer to find the best set 
of weights. Start point at (1/v,..) using bounds to limit the weights to 0..1. Then, 
it computes KL between the different probabilities. The overall complexity of opinion 
pooling finding is 
:math:`O[n^2 v^2]` +  :math:`O[w n^2]`.


2. **tSNE application**: this stage applies tSNE (t-distributed Stochastic Neighbor Embedding)
to the multiple views of the same data. Its complexity is about :math:`O[n^2 v m]`.

The overall complexity of MvtSNE is
:math:`O[n^2 v^2] + O[w n^2] + O[n^2 v m]`.

* :math:`n`: number of samples of each view
* :math:`v`: number of different views
* :math:`w`: weights dimension
* :math:`m`: input maximum iteration 

.. _mvsc:

MVSC
====

Given a multiview dataset with ``v`` input data matrices, multiview spectral clustering (MVSC) methods
produce a single clustering assignment, considering the information from all the 
input views.
Package ``multiview`` offers the class :class:`mvsc` to perform multiview spectral 
clustering. It is an extension to spectral clustering (``kernlab::specc`` in R language) 
to multiview datasets.

Complexity
----------
Multiview spectral clustering computes the diagonal matrix of the similarity matrices.
Firstly, computes distance matrices if necessary. After that, calculates laplacian matrices,
extracts commom principal components and apply KMeans algorithm to obtained data.
Roughly, the complexity of MVSC is
:math:`O[k n^2 v] + O[n \log(n)]`.

* :math:`n`: number of samples of each view
* :math:`k`: components to extract


Alternative use
===============

Although the methods in this package have been divided in dimensionality reduction
and clustering, there is a close relationship between both tasks. In fact, all three
methods can be used for both tasks.

First, the data projection produced by dimensionality reduction methods can be
fed to a standard clustering algorithm in order to obtain a multiview clustering.
Second, as mvsc also returns the projection resulting from the k first common
eigenvectors in matrix $evectors, this space can also be used as a low-dimensional
embedding of the original multiview data, for visualization or other purposes.

.. topic:: References: 

  * Abbas, Ali E. 2009. "A Kullback-Leibler View of Linear and Log-Linear Pools." *Decision Analysis*
    6 (1): 25-37. doi:`10.1287/deca.1080.0133 <http://pubsonline.informs.org/doi/abs/10.1287/deca.1080.0133>`_.

  * Carvalho, Arthur, and Kate Larson. 2012. "A Consensual Linear Opinion Pool." 
    http://arxiv.org/abs/1204.5399.

  * Kruskal, J B. 1964. "Multidimensional scaling by optimizing goodness of fit to 
    a nonmetric hypothesis." *Psychometrika* 29 (1): 1\-27. doi:`10.1007/BF02289565 <https://link.springer.com/article/10.1007%2FBF02289565>`_.

  * Ng, Andrew Y, Michael I Jordan, and Yair Weiss. 2001. "On spectral clustering: 
    Analysis and an algorithm." *Nips* 14 (14). MIT Press: 849-56.

  * Planck, Max, and Ulrike Von Luxburg. 2006. "A Tutorial on Spectral Clustering." 
    *Statistics and Computing* 17 (March). Springer US: 395-416. doi
    `10.1007/s11222-007-9033-z <https://link.springer.com/article/10.1007%2Fs11222-007-9033-z>`_.

  * Shi, Jianbo, and Jitendra Malik. 2005. "Normalized Cuts and Image Segmentation 
    Normalized Cuts and Image Segmentation." *Pattern Analysis and Machine Intelligence, IEEE Transactions* 
    on 22 (March): 888-905. doi:`10.1109/CVPR.1997.609407 <http://ieeexplore.ieee.org/document/609407/?reload=true>`_.

  * Trendafilov, Nickolay T. 2010. "Stepwise estimation of common principal 
    components." *Computational Statistics and Data Analysis* 54 (12): 3446-57. 
    doi:`10.1016/j.csda.2010.03.010 <http://www.sciencedirect.com/science/article/pii/S016794731000112X?via%3Dihub>`_.

  * Van Der Maaten, Laurens, Geoffrey Hinton, and Geoffrey Hinton van der Maaten. 
    2008. "Visualizing Data using t-SNE." doi:`10.1007/s10479-011-0841-3 <https://link.springer.com/article/10.1007%2Fs10479-011-0841-3>`_.

  * Multiview features dataset. https://archive.ics.uci.edu/ml/datasets/Multiple+Features