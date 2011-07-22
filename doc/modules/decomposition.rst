
=================================================================
Decomposing signals in components (matrix factorization problems)
=================================================================

.. currentmodule:: scikits.learn.decomposition


.. _PCA:

Principal component analysis (PCA)
==================================

Exact PCA and probabilistic interpretation
------------------------------------------

PCA is used to decompose a multivariate dataset in a set of successive
orthogonal components that explain a maximum amount of the variance. In
the scikit-learn, :class:`PCA` is implemented as a `transformer` object
that learns n components in its `fit` method, and can be used on new data
to project it on these components.

The optional parameter `whiten=True` parameter make it possible to
project the data onto the singular space while scaling each component
to unit variance. This is often useful if the models down-stream make
strong assumptions on the isotropy of the signal: this is for example
the case for Support Vector Machines with the RBF kernel and the K-Means
clustering algorithm. However in that case the inverse transform is no
longer exact since some information is lost while forward transforming.

In addition, the :class:`ProbabilisticPCA` object provides a
probabilistic interpretation of the PCA that can give a likelihood of
data based on the amount of variance it explains. As such it implements a
`score` method that can be used in cross-validation.

Below is an example of the iris dataset, which is comprised of 4
features, projected on the 2 dimensions that explain most variance:

.. figure:: ../auto_examples/decomposition/images/plot_pca_vs_lda_1.png
    :target: ../auto_examples/decomposition/plot_pca_vs_lda.html
    :align: center
    :scale: 75%

.. topic:: Examples:

    * :ref:`example_decomposition_plot_pca_vs_lda.py`


Approximate PCA
---------------

Often we are interested in projecting the data onto a lower dimensional
space that preserves most of the variance by dropping the singular vector
of components associated with lower singular values.

For instance for face recognition, if we work with 64x64 gray level pixel
pictures the dimensionality of the data is 4096 and it is slow to train a
RBF Support Vector Machine on such wide data. Furthermore we know that
intrinsic dimensionality of the data is much lower than 4096 since all
faces pictures look alike. The samples lie on a manifold of much lower
dimension (say around 200 for instance). The PCA algorithm can be used
to linearly transform the data while both reducing the dimensionality
and preserve most of the explained variance at the same time.

The class :class:`RandomizedPCA` is very useful in that case: since we
are going to drop most of the singular vectors it is much more efficient
to limit the computation to an approximated estimate of the singular
vectors we will keep to actually perform the transform.

:class:`RandomizedPCA` can hence be used as a drop in replacement for
:class:`PCA` minor the exception that we need to give it the size of
the lower dimensional space `n_components` as mandatory input parameter.

If we note :math:`n_{max} = max(n_{samples}, n_{features})` and
:math:`n_{min} = min(n_{samples}, n_{features})`, the time complexity
of :class:`RandomizedPCA` is :math:`O(n_{max}^2 \cdot n_{components})`
instead of :math:`O(n_{max}^2 \cdot n_{min})` for the exact method
implemented in :class:`PCA`.

The memory footprint of :class:`RandomizedPCA` is also proportional to
:math:`2 \cdot n_{max} \cdot n_{components}` instead of :math:`n_{max}
\cdot n_{min}` for the exact method.

Furthermore :class:`RandomizedPCA` is able to work with
`scipy.sparse` matrices as input which make it suitable for reducing
the dimensionality of features extracted from text documents for
instance.

Note: the implementation of `inverse_transform` in :class:`RandomizedPCA`
is not the exact inverse transform of `transform` even when
`whiten=False` (default).


.. topic:: Examples:

    * :ref:`example_applications_face_recognition.py`
    * :ref:`example_decomposition_plot_faces_decomposition.py`

.. topic:: References:

    * `"Finding structure with randomness: Stochastic algorithms for
      constructing approximate matrix decompositions"
      <http://arxiv.org/abs/0909.4061>`_
      Halko, et al., 2009

.. _kernel_PCA:

Kernel PCA
----------

:class:`KernelPCA` is an extension of PCA which achieves non-linear
dimensionality reduction through the use of kernels. It has many
applications including denoising, compression and structured prediction
(kernel dependency estimation). :class:`KernelPCA` supports both
`transform` and `inverse_transform`.

.. figure:: ../auto_examples/decomposition/images/plot_kernel_pca_1.png
    :target: ../auto_examples/decomposition/plot_kernel_pca.html
    :align: center
    :scale: 75%

.. topic:: Examples:

    * :ref:`example_decomposition_plot_kernel_pca.py`

.. _SparsePCA:

Sparse Principal Components Analysis (SparsePCA)
------------------------------------------------

:class:`SparsePCA` is a variant of PCA, with the goal of extracting the
set of sparse components that best reconstruct the data.

Principal component analysis (:class:`PCA`) has the disadvantage that the
components extracted by this method have exclusively dense expressions, i.e.
they have non-zero coefficients when expressed as linear combinations of the
original variables. This can make interpretation difficult. In many cases,
the real underlying components can be more naturally imagined as sparse
vectors; for example in face recognition, components might naturally map to
parts of faces.

Sparse principal components yields a more parsimonious, interpretable
representation, clearly emphasizing which of the original features contribute
to the differences between samples.

The following example illustrates 12 components extracted using sparse PCA
with a value of `alpha=5` on the digits dataset. Only images of the digit 3
were considered.  It can be seen how the regularization term induces many
zeros. Furthermore, the natural structure of the data causes the non-zero
coefficients to be vertically adjacent. The model does not enforce this
mathematically: each component is a vector :math:`h \in \mathbf{R}^{64}`, and there
is no notion of vertical adjacency except during the human-friendly
visualization as 8x8 pixel images.
The fact that the components shown below appear local
is the effect of the inherent structure of the data, which makes such local
patterns minimize reconstruction error. There exist sparsity-inducing norms
that take into account adjacency and different kinds of structure; see
see [Jen09] for a review of such methods. For more details on how to use
Sparse PCA, see the `Examples` section below.


.. figure:: ../auto_examples/decomposition/images/plot_faces_decomposition_4.png
   :target: ../auto_examples/decomposition/plot_faces_decomposition.html
   :align: center
   :scale: 50%


Note that there are many different formulations for the Sparse PCA
problem. The one implemented here is based on [Mrl09]_ . The optimization
problem solved is a PCA problem (dictionary learning) with an
:math:`\ell_1` penalty on the components:

.. math::
   (U^*, V^*) = \underset{U, V}{\operatorname{arg\,min\,}} & \frac{1}{2}
                ||X-UV||_2^2+\alpha||V||_1 \\
                \text{subject to\,} & ||U_k||_2 = 1 \text{ for all }
                0 \leq k < n_{components}


The sparsity inducing :math:`\ell_1` norm also prevents learning
components from noise when few training samples are available. The degree
of penalization (and thus sparsity) can be adjusted through the
hyperparameter `alpha`. Small values lead to a gently regularized
factorization, while larger values shrink many coefficients to zero.


.. topic:: Examples:

   * :ref:`example_decomposition_plot_faces_decomposition.py`

.. topic:: References:

   * [Mrl09] `"Online Dictionary Learning for Sparse Coding"
     <http://www.di.ens.fr/sierra/pdfs/icml09.pdf>`_
     J. Mairal, F. Bach, J. Ponce, G. Sapiro, 2009
   * [Jen09] `"Structured Sparse Principal Component Analysis"
     <www.di.ens.fr/~fbach/sspca_AISTATS2010.pdf>`_
     R. Jenatton, G. Obozinski, F. Bach, 2009


.. _MiniBatchSparsePCA:

Mini Batch Sparse Principal Components Analysis (MiniBatchSparsePCA)
--------------------------------------------------------------------

Mini Batch Sparse PCA (:class:`MiniBatchSparsePCA`) is a variant of
:class:`SparsePCA` that is faster but less accurate. The increased speed is
reached by iterating over small chunks of the set of features, for a given
number of iterations.

Note that while this is in the spirit of an online algorithm, the class
:class:`MiniBatchSparsePCA` does not implement `partial_fit` because the
algorithm is online along the features direction, not the samples direction.


.. topic:: Examples:

   * :ref:`example_decomposition_plot_faces_decomposition.py`


.. topic:: References:

   * [Mrl09] `"Online Dictionary Learning for Sparse Coding"
     <http://www.di.ens.fr/sierra/pdfs/icml09.pdf>`_
     J. Mairal, F. Bach, J. Ponce, G. Sapiro, 2009
.. _ICA:

Independent component analysis (ICA)
====================================

Independent component analysis separates a multivariate signal into
additive subcomponents that are maximally independent. It is
implemented in scikit-learn using the :class:`Fast ICA <FastICA>`
algorithm.

It is classically used to separate mixed signals (a problem known as
*blind source separation*), as in the example below:

.. figure:: ../auto_examples/decomposition/images/plot_ica_blind_source_separation_1.png
    :target: ../auto_examples/decomposition/plot_ica_blind_source_separation.html
    :align: center
    :scale: 50%


.. topic:: Examples:

    * :ref:`example_decomposition_plot_ica_blind_source_separation.py`
    * :ref:`example_decomposition_plot_ica_vs_pca.py`
    * :ref:`example_decomposition_plot_digits_decomposition.py`


.. _NMF:

Non-negative matrix factorization (NMF)
=======================================

:class:`NMF` is an alternative approach to decomposition that assumes that the
data and the components are non-negative. :class:`NMF` can be plugged in
instead of :class:`PCA` or its variants, in the cases where the data matrix
does not contain negative values.

Unlike :class:`PCA`, the representation of a vector is obtained in an additive
fashion, by superimposing the components, without substracting. Such additive
models are efficient for representing images and text.

It has been observed in [Hoyer, 04] that, when carefully constrained,
:class:`NMF` can produce a parts-based representation of the dataset,
resulting in interpretable models. The following example displays 16
sparse components found by :class:`NMF` on the images of the digit 3 from the
digits dataset.

.. |pca_img| image:: ../auto_examples/decomposition/images/plot_faces_decomposition_1.png
    :target: ../auto_examples/decomposition/plot_faces_decomposition.html
    :scale: 50%

.. |nmf_img| image:: ../auto_examples/decomposition/images/plot_faces_decomposition_2.png
    :target: ../auto_examples/decomposition/plot_faces_decomposition.html
    :scale: 50%

.. centered:: |pca_img| |nmf_img|


The :attr:`init` attribute determines the initialization method applied, which
has a great impact on the performance of the method. :class:`NMF` implements
the method Nonnegative Double Singular Value Decomposition. NNDSVD is based on
two SVD processes, one approximating the data matrix, the other approximating
positive sections of the resulting partial SVD factors utilizing an algebraic
property of unit rank matrices. The basic NNDSVD algorithm is better fit for
sparse factorization. Its variants NNDSVDa (in which all zeros are set equal to
the mean of all elements of the data), and NNDSVDar (in which the zeros are set
to random perturbations less than the mean of the data divided by 100) are
recommended in the dense case.

:class:`NMF` can also be initialized with random non-negative matrices, by
passing an integer seed or a `RandomState` to :attr:`init`.

In :class:`NMF`, sparseness can be enforced by setting the attribute
:attr:`sparseness` to `data` or `components`. Sparse components lead to
localized features, and sparse data leads to a more efficient representation
of the data.

.. topic:: Examples:

    * :ref:`example_decomposition_plot_faces_decomposition.py`

.. topic:: References:

    * `"Learning the parts of objects by non-negative matrix factorization"
      <http://www.seas.upenn.edu/~ddlee/Papers/nmf.pdf>`_
      D. Lee, S. Seung, 1999

    * `"Non-negative Matrix Factorization with Sparseness Constraints"
      <http://www.cs.helsinki.fi/u/phoyer/papers/pdf/NMFscweb.pdf>`_
      P. Hoyer, 2004

    * `"Projected gradient methods for non-negative matrix factorization"
      <http://www.csie.ntu.edu.tw/~cjlin/nmf/>`_
      C.-J. Lin, 2007

    * `"SVD based initialization: A head start for nonnegative
      matrix factorization"
      <http://www.cs.rpi.edu/~boutsc/files/nndsvd.pdf>`_
      C. Boutsidis, E. Gallopoulos, 2008
