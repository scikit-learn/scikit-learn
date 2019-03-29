.. _decompositions:


=================================================================
Decomposing signals in components (matrix factorization problems)
=================================================================

.. currentmodule:: sklearn.decomposition


.. _PCA:


Principal component analysis (PCA)
==================================

Exact PCA and probabilistic interpretation
------------------------------------------

PCA is used to decompose a multivariate dataset in a set of successive
orthogonal components that explain a maximum amount of the variance. In
scikit-learn, :class:`PCA` is implemented as a *transformer* object
that learns :math:`n` components in its ``fit`` method, and can be used on new
data to project it on these components.

PCA centers but does not scale the input data for each feature before
applying the SVD. The optional parameter parameter ``whiten=True`` makes it
possible to project the data onto the singular space while scaling each
component to unit variance. This is often useful if the models down-stream make
strong assumptions on the isotropy of the signal: this is for example the case
for Support Vector Machines with the RBF kernel and the K-Means clustering
algorithm.

Below is an example of the iris dataset, which is comprised of 4
features, projected on the 2 dimensions that explain most variance:

.. figure:: ../auto_examples/decomposition/images/sphx_glr_plot_pca_vs_lda_001.png
    :target: ../auto_examples/decomposition/plot_pca_vs_lda.html
    :align: center
    :scale: 75%


The :class:`PCA` object also provides a
probabilistic interpretation of the PCA that can give a likelihood of
data based on the amount of variance it explains. As such it implements a
`score` method that can be used in cross-validation:

.. figure:: ../auto_examples/decomposition/images/sphx_glr_plot_pca_vs_fa_model_selection_001.png
    :target: ../auto_examples/decomposition/plot_pca_vs_fa_model_selection.html
    :align: center
    :scale: 75%


.. topic:: Examples:

    * :ref:`sphx_glr_auto_examples_decomposition_plot_pca_vs_lda.py`
    * :ref:`sphx_glr_auto_examples_decomposition_plot_pca_vs_fa_model_selection.py`


.. _IncrementalPCA:

Incremental PCA
---------------

The :class:`PCA` object is very useful, but has certain limitations for
large datasets. The biggest limitation is that :class:`PCA` only supports
batch processing, which means all of the data to be processed must fit in main
memory. The :class:`IncrementalPCA` object uses a different form of
processing and allows for partial computations which almost
exactly match the results of :class:`PCA` while processing the data in a
minibatch fashion. :class:`IncrementalPCA` makes it possible to implement
out-of-core Principal Component Analysis either by:

 * Using its ``partial_fit`` method on chunks of data fetched sequentially
   from the local hard drive or a network database.

 * Calling its fit method on a memory mapped file using ``numpy.memmap``.

:class:`IncrementalPCA` only stores estimates of component and noise variances,
in order update ``explained_variance_ratio_`` incrementally. This is why
memory usage depends on the number of samples per batch, rather than the
number of samples to be processed in the dataset.

As in :class:`PCA`, :class:`IncrementalPCA` centers but does not scale the
input data for each feature before applying the SVD.

.. figure:: ../auto_examples/decomposition/images/sphx_glr_plot_incremental_pca_001.png
    :target: ../auto_examples/decomposition/plot_incremental_pca.html
    :align: center
    :scale: 75%

.. figure:: ../auto_examples/decomposition/images/sphx_glr_plot_incremental_pca_002.png
    :target: ../auto_examples/decomposition/plot_incremental_pca.html
    :align: center
    :scale: 75%


.. topic:: Examples:

    * :ref:`sphx_glr_auto_examples_decomposition_plot_incremental_pca.py`


.. _RandomizedPCA:

PCA using randomized SVD
------------------------

It is often interesting to project data to a lower-dimensional
space that preserves most of the variance, by dropping the singular vector
of components associated with lower singular values.

For instance, if we work with 64x64 pixel gray-level pictures
for face recognition,
the dimensionality of the data is 4096 and it is slow to train an
RBF support vector machine on such wide data. Furthermore we know that
the intrinsic dimensionality of the data is much lower than 4096 since all
pictures of human faces look somewhat alike.
The samples lie on a manifold of much lower
dimension (say around 200 for instance). The PCA algorithm can be used
to linearly transform the data while both reducing the dimensionality
and preserve most of the explained variance at the same time.

The class :class:`PCA` used with the optional parameter
``svd_solver='randomized'`` is very useful in that case: since we are going
to drop most of the singular vectors it is much more efficient to limit the
computation to an approximated estimate of the singular vectors we will keep
to actually perform the transform.

For instance, the following shows 16 sample portraits (centered around
0.0) from the Olivetti dataset. On the right hand side are the first 16
singular vectors reshaped as portraits. Since we only require the top
16 singular vectors of a dataset with size :math:`n_{samples} = 400`
and :math:`n_{features} = 64 \times 64 = 4096`, the computation time is
less than 1s:

.. |orig_img| image:: ../auto_examples/decomposition/images/sphx_glr_plot_faces_decomposition_001.png
   :target: ../auto_examples/decomposition/plot_faces_decomposition.html
   :scale: 60%

.. |pca_img| image:: ../auto_examples/decomposition/images/sphx_glr_plot_faces_decomposition_002.png
   :target: ../auto_examples/decomposition/plot_faces_decomposition.html
   :scale: 60%

.. centered:: |orig_img| |pca_img|

If we note :math:`n_{\max} = \max(n_{\mathrm{samples}}, n_{\mathrm{features}})` and
:math:`n_{\min} = \min(n_{\mathrm{samples}}, n_{\mathrm{features}})`, the time complexity
of the randomized :class:`PCA` is :math:`O(n_{\max}^2 \cdot n_{\mathrm{components}})`
instead of :math:`O(n_{\max}^2 \cdot n_{\min})` for the exact method
implemented in :class:`PCA`.

The memory footprint of randomized :class:`PCA` is also proportional to
:math:`2 \cdot n_{\max} \cdot n_{\mathrm{components}}` instead of :math:`n_{\max}
\cdot n_{\min}` for the exact method.

Note: the implementation of ``inverse_transform`` in :class:`PCA` with
``svd_solver='randomized'`` is not the exact inverse transform of
``transform`` even when ``whiten=False`` (default).


.. topic:: Examples:

    * :ref:`sphx_glr_auto_examples_applications_plot_face_recognition.py`
    * :ref:`sphx_glr_auto_examples_decomposition_plot_faces_decomposition.py`

.. topic:: References:

    * `"Finding structure with randomness: Stochastic algorithms for
      constructing approximate matrix decompositions"
      <https://arxiv.org/abs/0909.4061>`_
      Halko, et al., 2009


.. _kernel_PCA:

Kernel PCA
----------

:class:`KernelPCA` is an extension of PCA which achieves non-linear
dimensionality reduction through the use of kernels (see :ref:`metrics`). It
has many applications including denoising, compression and structured
prediction (kernel dependency estimation). :class:`KernelPCA` supports both
``transform`` and ``inverse_transform``.

.. figure:: ../auto_examples/decomposition/images/sphx_glr_plot_kernel_pca_001.png
    :target: ../auto_examples/decomposition/plot_kernel_pca.html
    :align: center
    :scale: 75%

.. topic:: Examples:

    * :ref:`sphx_glr_auto_examples_decomposition_plot_kernel_pca.py`


.. _SparsePCA:

Sparse principal components analysis (SparsePCA and MiniBatchSparsePCA)
-----------------------------------------------------------------------

:class:`SparsePCA` is a variant of PCA, with the goal of extracting the
set of sparse components that best reconstruct the data.

Mini-batch sparse PCA (:class:`MiniBatchSparsePCA`) is a variant of
:class:`SparsePCA` that is faster but less accurate. The increased speed is
reached by iterating over small chunks of the set of features, for a given
number of iterations.


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

The following example illustrates 16 components extracted using sparse PCA from
the Olivetti faces dataset.  It can be seen how the regularization term induces
many zeros. Furthermore, the natural structure of the data causes the non-zero
coefficients to be vertically adjacent. The model does not enforce this
mathematically: each component is a vector :math:`h \in \mathbf{R}^{4096}`, and
there is no notion of vertical adjacency except during the human-friendly
visualization as 64x64 pixel images. The fact that the components shown below
appear local is the effect of the inherent structure of the data, which makes
such local patterns minimize reconstruction error. There exist sparsity-inducing
norms that take into account adjacency and different kinds of structure; see
[Jen09]_ for a review of such methods.
For more details on how to use Sparse PCA, see the Examples section, below.


.. |spca_img| image:: ../auto_examples/decomposition/images/sphx_glr_plot_faces_decomposition_005.png
   :target: ../auto_examples/decomposition/plot_faces_decomposition.html
   :scale: 60%

.. centered:: |pca_img| |spca_img|

Note that there are many different formulations for the Sparse PCA
problem. The one implemented here is based on [Mrl09]_ . The optimization
problem solved is a PCA problem (dictionary learning) with an
:math:`\ell_1` penalty on the components:

.. math::
   (U^*, V^*) = \underset{U, V}{\operatorname{arg\,min\,}} & \frac{1}{2}
                ||X-UV||_2^2+\alpha||V||_1 \\
                \text{subject to } & ||U_k||_2 = 1 \text{ for all }
                0 \leq k < n_{components}


The sparsity-inducing :math:`\ell_1` norm also prevents learning
components from noise when few training samples are available. The degree
of penalization (and thus sparsity) can be adjusted through the
hyperparameter ``alpha``. Small values lead to a gently regularized
factorization, while larger values shrink many coefficients to zero.

.. note::

  While in the spirit of an online algorithm, the class
  :class:`MiniBatchSparsePCA` does not implement ``partial_fit`` because
  the algorithm is online along the features direction, not the samples
  direction.

.. topic:: Examples:

   * :ref:`sphx_glr_auto_examples_decomposition_plot_faces_decomposition.py`

.. topic:: References:

  .. [Mrl09] `"Online Dictionary Learning for Sparse Coding"
     <https://www.di.ens.fr/sierra/pdfs/icml09.pdf>`_
     J. Mairal, F. Bach, J. Ponce, G. Sapiro, 2009
  .. [Jen09] `"Structured Sparse Principal Component Analysis"
     <https://www.di.ens.fr/~fbach/sspca_AISTATS2010.pdf>`_
     R. Jenatton, G. Obozinski, F. Bach, 2009


.. _LSA:

Truncated singular value decomposition and latent semantic analysis
===================================================================

:class:`TruncatedSVD` implements a variant of singular value decomposition
(SVD) that only computes the :math:`k` largest singular values,
where :math:`k` is a user-specified parameter.

When truncated SVD is applied to term-document matrices
(as returned by ``CountVectorizer`` or ``TfidfVectorizer``),
this transformation is known as
`latent semantic analysis <https://nlp.stanford.edu/IR-book/pdf/18lsi.pdf>`_
(LSA), because it transforms such matrices
to a "semantic" space of low dimensionality.
In particular, LSA is known to combat the effects of synonymy and polysemy
(both of which roughly mean there are multiple meanings per word),
which cause term-document matrices to be overly sparse
and exhibit poor similarity under measures such as cosine similarity.

.. note::
    LSA is also known as latent semantic indexing, LSI,
    though strictly that refers to its use in persistent indexes
    for information retrieval purposes.

Mathematically, truncated SVD applied to training samples :math:`X`
produces a low-rank approximation :math:`X`:

.. math::
    X \approx X_k = U_k \Sigma_k V_k^\top

After this operation, :math:`U_k \Sigma_k^\top`
is the transformed training set with :math:`k` features
(called ``n_components`` in the API).

To also transform a test set :math:`X`, we multiply it with :math:`V_k`:

.. math::
    X' = X V_k

.. note::
    Most treatments of LSA in the natural language processing (NLP)
    and information retrieval (IR) literature
    swap the axes of the matrix :math:`X` so that it has shape
    ``n_features`` × ``n_samples``.
    We present LSA in a different way that matches the scikit-learn API better,
    but the singular values found are the same.

:class:`TruncatedSVD` is very similar to :class:`PCA`, but differs
in that it works on sample matrices :math:`X` directly
instead of their covariance matrices.
When the columnwise (per-feature) means of :math:`X`
are subtracted from the feature values,
truncated SVD on the resulting matrix is equivalent to PCA.
In practical terms, this means
that the :class:`TruncatedSVD` transformer accepts ``scipy.sparse``
matrices without the need to densify them,
as densifying may fill up memory even for medium-sized document collections.

While the :class:`TruncatedSVD` transformer
works with any (sparse) feature matrix,
using it on tf–idf matrices is recommended over raw frequency counts
in an LSA/document processing setting.
In particular, sublinear scaling and inverse document frequency
should be turned on (``sublinear_tf=True, use_idf=True``)
to bring the feature values closer to a Gaussian distribution,
compensating for LSA's erroneous assumptions about textual data.

.. topic:: Examples:

   * :ref:`sphx_glr_auto_examples_text_plot_document_clustering.py`

.. topic:: References:

  * Christopher D. Manning, Prabhakar Raghavan and Hinrich Schütze (2008),
    *Introduction to Information Retrieval*, Cambridge University Press,
    chapter 18: `Matrix decompositions & latent semantic indexing
    <https://nlp.stanford.edu/IR-book/pdf/18lsi.pdf>`_


.. _DictionaryLearning:

Dictionary Learning
===================

.. _SparseCoder:

Sparse coding with a precomputed dictionary
-------------------------------------------

The :class:`SparseCoder` object is an estimator that can be used to transform signals
into sparse linear combination of atoms from a fixed, precomputed dictionary
such as a discrete wavelet basis. This object therefore does not
implement a ``fit`` method. The transformation amounts
to a sparse coding problem: finding a representation of the data as a linear
combination of as few dictionary atoms as possible. All variations of
dictionary learning implement the following transform methods, controllable via
the ``transform_method`` initialization parameter:

* Orthogonal matching pursuit (:ref:`omp`)

* Least-angle regression (:ref:`least_angle_regression`)

* Lasso computed by least-angle regression

* Lasso using coordinate descent (:ref:`lasso`)

* Thresholding

Thresholding is very fast but it does not yield accurate reconstructions.
They have been shown useful in literature for classification tasks. For image
reconstruction tasks, orthogonal matching pursuit yields the most accurate,
unbiased reconstruction.

The dictionary learning objects offer, via the ``split_code`` parameter, the
possibility to separate the positive and negative values in the results of
sparse coding. This is useful when dictionary learning is used for extracting
features that will be used for supervised learning, because it allows the
learning algorithm to assign different weights to negative loadings of a
particular atom, from to the corresponding positive loading.

The split code for a single sample has length ``2 * n_components``
and is constructed using the following rule: First, the regular code of length
``n_components`` is computed. Then, the first ``n_components`` entries of the
``split_code`` are
filled with the positive part of the regular code vector. The second half of
the split code is filled with the negative part of the code vector, only with
a positive sign. Therefore, the split_code is non-negative.


.. topic:: Examples:

    * :ref:`sphx_glr_auto_examples_decomposition_plot_sparse_coding.py`


Generic dictionary learning
---------------------------

Dictionary learning (:class:`DictionaryLearning`) is a matrix factorization
problem that amounts to finding a (usually overcomplete) dictionary that will
perform well at sparsely encoding the fitted data.

Representing data as sparse combinations of atoms from an overcomplete
dictionary is suggested to be the way the mammalian primary visual cortex works.
Consequently, dictionary learning applied on image patches has been shown to
give good results in image processing tasks such as image completion,
inpainting and denoising, as well as for supervised recognition tasks.

Dictionary learning is an optimization problem solved by alternatively updating
the sparse code, as a solution to multiple Lasso problems, considering the
dictionary fixed, and then updating the dictionary to best fit the sparse code.

.. math::
   (U^*, V^*) = \underset{U, V}{\operatorname{arg\,min\,}} & \frac{1}{2}
                ||X-UV||_2^2+\alpha||U||_1 \\
                \text{subject to } & ||V_k||_2 = 1 \text{ for all }
                0 \leq k < n_{\mathrm{atoms}}


.. |pca_img2| image:: ../auto_examples/decomposition/images/sphx_glr_plot_faces_decomposition_002.png
   :target: ../auto_examples/decomposition/plot_faces_decomposition.html
   :scale: 60%

.. |dict_img2| image:: ../auto_examples/decomposition/images/sphx_glr_plot_faces_decomposition_006.png
   :target: ../auto_examples/decomposition/plot_faces_decomposition.html
   :scale: 60%

.. centered:: |pca_img2| |dict_img2|


After using such a procedure to fit the dictionary, the transform is simply a
sparse coding step that shares the same implementation with all dictionary
learning objects (see :ref:`SparseCoder`).

It is also possible to constrain the dictionary and/or code to be positive to
match constraints that may be present in the data. Below are the faces with
different positivity constraints applied. Red indicates negative values, blue
indicates positive values, and white represents zeros.


.. |dict_img_pos1| image:: ../auto_examples/decomposition/images/sphx_glr_plot_faces_decomposition_011.png
    :target: ../auto_examples/decomposition/plot_image_denoising.html
    :scale: 60%

.. |dict_img_pos2| image:: ../auto_examples/decomposition/images/sphx_glr_plot_faces_decomposition_012.png
    :target: ../auto_examples/decomposition/plot_image_denoising.html
    :scale: 60%

.. |dict_img_pos3| image:: ../auto_examples/decomposition/images/sphx_glr_plot_faces_decomposition_013.png
    :target: ../auto_examples/decomposition/plot_image_denoising.html
    :scale: 60%

.. |dict_img_pos4| image:: ../auto_examples/decomposition/images/sphx_glr_plot_faces_decomposition_014.png
    :target: ../auto_examples/decomposition/plot_image_denoising.html
    :scale: 60%

.. centered:: |dict_img_pos1| |dict_img_pos2|
.. centered:: |dict_img_pos3| |dict_img_pos4|


The following image shows how a dictionary learned from 4x4 pixel image patches
extracted from part of the image of a raccoon face looks like.


.. figure:: ../auto_examples/decomposition/images/sphx_glr_plot_image_denoising_001.png
    :target: ../auto_examples/decomposition/plot_image_denoising.html
    :align: center
    :scale: 50%


.. topic:: Examples:

  * :ref:`sphx_glr_auto_examples_decomposition_plot_image_denoising.py`


.. topic:: References:

  * `"Online dictionary learning for sparse coding"
    <https://www.di.ens.fr/sierra/pdfs/icml09.pdf>`_
    J. Mairal, F. Bach, J. Ponce, G. Sapiro, 2009

.. _MiniBatchDictionaryLearning:

Mini-batch dictionary learning
------------------------------

:class:`MiniBatchDictionaryLearning` implements a faster, but less accurate
version of the dictionary learning algorithm that is better suited for large
datasets.

By default, :class:`MiniBatchDictionaryLearning` divides the data into
mini-batches and optimizes in an online manner by cycling over the mini-batches
for the specified number of iterations. However, at the moment it does not
implement a stopping condition.

The estimator also implements ``partial_fit``, which updates the dictionary by
iterating only once over a mini-batch. This can be used for online learning
when the data is not readily available from the start, or for when the data
does not fit into the memory.

.. currentmodule:: sklearn.cluster

.. image:: ../auto_examples/cluster/images/sphx_glr_plot_dict_face_patches_001.png
    :target: ../auto_examples/cluster/plot_dict_face_patches.html
    :scale: 50%
    :align: right

.. topic:: **Clustering for dictionary learning**

   Note that when using dictionary learning to extract a representation
   (e.g. for sparse coding) clustering can be a good proxy to learn the
   dictionary. For instance the :class:`MiniBatchKMeans` estimator is
   computationally efficient and implements on-line learning with a
   ``partial_fit`` method.

    Example: :ref:`sphx_glr_auto_examples_cluster_plot_dict_face_patches.py`

.. currentmodule:: sklearn.decomposition

.. _FA:

Factor Analysis
===============

In unsupervised learning we only have a dataset :math:`X = \{x_1, x_2, \dots, x_n
\}`. How can this dataset be described mathematically? A very simple
`continuous latent variable` model for :math:`X` is

.. math:: x_i = W h_i + \mu + \epsilon

The vector :math:`h_i` is called "latent" because it is unobserved. :math:`\epsilon` is
considered a noise term distributed according to a Gaussian with mean 0 and
covariance :math:`\Psi` (i.e. :math:`\epsilon \sim \mathcal{N}(0, \Psi)`), :math:`\mu` is some
arbitrary offset vector. Such a model is called "generative" as it describes
how :math:`x_i` is generated from :math:`h_i`. If we use all the :math:`x_i`'s as columns to form
a matrix :math:`\mathbf{X}` and all the :math:`h_i`'s as columns of a matrix :math:`\mathbf{H}`
then we can write (with suitably defined :math:`\mathbf{M}` and :math:`\mathbf{E}`):

.. math::
    \mathbf{X} = W \mathbf{H} + \mathbf{M} + \mathbf{E}

In other words, we *decomposed* matrix :math:`\mathbf{X}`.

If :math:`h_i` is given, the above equation automatically implies the following
probabilistic interpretation:

.. math:: p(x_i|h_i) = \mathcal{N}(Wh_i + \mu, \Psi)

For a complete probabilistic model we also need a prior distribution for the
latent variable :math:`h`. The most straightforward assumption (based on the nice
properties of the Gaussian distribution) is :math:`h \sim \mathcal{N}(0,
\mathbf{I})`.  This yields a Gaussian as the marginal distribution of :math:`x`:

.. math:: p(x) = \mathcal{N}(\mu, WW^T + \Psi)

Now, without any further assumptions the idea of having a latent variable :math:`h`
would be superfluous -- :math:`x` can be completely modelled with a mean
and a covariance. We need to impose some more specific structure on one
of these two parameters. A simple additional assumption regards the
structure of the error covariance :math:`\Psi`:

* :math:`\Psi = \sigma^2 \mathbf{I}`: This assumption leads to
  the probabilistic model of :class:`PCA`.

* :math:`\Psi = \mathrm{diag}(\psi_1, \psi_2, \dots, \psi_n)`: This model is called
  :class:`FactorAnalysis`, a classical statistical model. The matrix W is
  sometimes called the "factor loading matrix".

Both models essentially estimate a Gaussian with a low-rank covariance matrix.
Because both models are probabilistic they can be integrated in more complex
models, e.g. Mixture of Factor Analysers. One gets very different models (e.g.
:class:`FastICA`) if non-Gaussian priors on the latent variables are assumed.

Factor analysis *can* produce similar components (the columns of its loading
matrix) to :class:`PCA`. However, one can not make any general statements
about these components (e.g. whether they are orthogonal):

.. |pca_img3| image:: ../auto_examples/decomposition/images/sphx_glr_plot_faces_decomposition_002.png
    :target: ../auto_examples/decomposition/plot_faces_decomposition.html
    :scale: 60%

.. |fa_img3| image:: ../auto_examples/decomposition/images/sphx_glr_plot_faces_decomposition_009.png
    :target: ../auto_examples/decomposition/plot_faces_decomposition.html
    :scale: 60%

.. centered:: |pca_img3| |fa_img3|

The main advantage for Factor Analysis over :class:`PCA` is that
it can model the variance in every direction of the input space independently
(heteroscedastic noise):

.. figure:: ../auto_examples/decomposition/images/sphx_glr_plot_faces_decomposition_008.png
    :target: ../auto_examples/decomposition/plot_faces_decomposition.html
    :align: center
    :scale: 75%

This allows better model selection than probabilistic PCA in the presence
of heteroscedastic noise:

.. figure:: ../auto_examples/decomposition/images/sphx_glr_plot_pca_vs_fa_model_selection_002.png
    :target: ../auto_examples/decomposition/plot_pca_vs_fa_model_selection.html
    :align: center
    :scale: 75%


.. topic:: Examples:

    * :ref:`sphx_glr_auto_examples_decomposition_plot_pca_vs_fa_model_selection.py`

.. _ICA:

Independent component analysis (ICA)
====================================

Independent component analysis separates a multivariate signal into
additive subcomponents that are maximally independent. It is
implemented in scikit-learn using the :class:`Fast ICA <FastICA>`
algorithm. Typically, ICA is not used for reducing dimensionality but
for separating superimposed signals. Since the ICA model does not include
a noise term, for the model to be correct, whitening must be applied.
This can be done internally using the whiten argument or manually using one
of the PCA variants.

It is classically used to separate mixed signals (a problem known as
*blind source separation*), as in the example below:

.. figure:: ../auto_examples/decomposition/images/sphx_glr_plot_ica_blind_source_separation_001.png
    :target: ../auto_examples/decomposition/plot_ica_blind_source_separation.html
    :align: center
    :scale: 60%


ICA can also be used as yet another non linear decomposition that finds
components with some sparsity:

.. |pca_img4| image:: ../auto_examples/decomposition/images/sphx_glr_plot_faces_decomposition_002.png
    :target: ../auto_examples/decomposition/plot_faces_decomposition.html
    :scale: 60%

.. |ica_img4| image:: ../auto_examples/decomposition/images/sphx_glr_plot_faces_decomposition_004.png
    :target: ../auto_examples/decomposition/plot_faces_decomposition.html
    :scale: 60%

.. centered:: |pca_img4| |ica_img4|

.. topic:: Examples:

    * :ref:`sphx_glr_auto_examples_decomposition_plot_ica_blind_source_separation.py`
    * :ref:`sphx_glr_auto_examples_decomposition_plot_ica_vs_pca.py`
    * :ref:`sphx_glr_auto_examples_decomposition_plot_faces_decomposition.py`


.. _NMF:

Non-negative matrix factorization (NMF or NNMF)
===============================================

NMF with the Frobenius norm
---------------------------

:class:`NMF` [1]_ is an alternative approach to decomposition that assumes that the
data and the components are non-negative. :class:`NMF` can be plugged in
instead of :class:`PCA` or its variants, in the cases where the data matrix
does not contain negative values. It finds a decomposition of samples
:math:`X` into two matrices :math:`W` and :math:`H` of non-negative elements,
by optimizing the distance :math:`d` between :math:`X` and the matrix product
:math:`WH`. The most widely used distance function is the squared Frobenius
norm, which is an obvious extension of the Euclidean norm to matrices:

.. math::
    d_{\mathrm{Fro}}(X, Y) = \frac{1}{2} ||X - Y||_{\mathrm{Fro}}^2 = \frac{1}{2} \sum_{i,j} (X_{ij} - {Y}_{ij})^2

Unlike :class:`PCA`, the representation of a vector is obtained in an additive
fashion, by superimposing the components, without subtracting. Such additive
models are efficient for representing images and text.

It has been observed in [Hoyer, 2004] [2]_ that, when carefully constrained,
:class:`NMF` can produce a parts-based representation of the dataset,
resulting in interpretable models. The following example displays 16
sparse components found by :class:`NMF` from the images in the Olivetti
faces dataset, in comparison with the PCA eigenfaces.

.. |pca_img5| image:: ../auto_examples/decomposition/images/sphx_glr_plot_faces_decomposition_002.png
    :target: ../auto_examples/decomposition/plot_faces_decomposition.html
    :scale: 60%

.. |nmf_img5| image:: ../auto_examples/decomposition/images/sphx_glr_plot_faces_decomposition_003.png
    :target: ../auto_examples/decomposition/plot_faces_decomposition.html
    :scale: 60%

.. centered:: |pca_img5| |nmf_img5|


The :attr:`init` attribute determines the initialization method applied, which
has a great impact on the performance of the method. :class:`NMF` implements the
method Nonnegative Double Singular Value Decomposition. NNDSVD [4]_ is based on
two SVD processes, one approximating the data matrix, the other approximating
positive sections of the resulting partial SVD factors utilizing an algebraic
property of unit rank matrices. The basic NNDSVD algorithm is better fit for
sparse factorization. Its variants NNDSVDa (in which all zeros are set equal to
the mean of all elements of the data), and NNDSVDar (in which the zeros are set
to random perturbations less than the mean of the data divided by 100) are
recommended in the dense case.

Note that the Multiplicative Update ('mu') solver cannot update zeros present in
the initialization, so it leads to poorer results when used jointly with the
basic NNDSVD algorithm which introduces a lot of zeros; in this case, NNDSVDa or
NNDSVDar should be preferred.

:class:`NMF` can also be initialized with correctly scaled random non-negative
matrices by setting :attr:`init="random"`. An integer seed or a
``RandomState`` can also be passed to :attr:`random_state` to control
reproducibility.

In :class:`NMF`, L1 and L2 priors can be added to the loss function in order
to regularize the model. The L2 prior uses the Frobenius norm, while the L1
prior uses an elementwise L1 norm. As in :class:`ElasticNet`, we control the
combination of L1 and L2 with the :attr:`l1_ratio` (:math:`\rho`) parameter,
and the intensity of the regularization with the :attr:`alpha`
(:math:`\alpha`) parameter. Then the priors terms are:

.. math::
    \alpha \rho ||W||_1 + \alpha \rho ||H||_1
    + \frac{\alpha(1-\rho)}{2} ||W||_{\mathrm{Fro}} ^ 2
    + \frac{\alpha(1-\rho)}{2} ||H||_{\mathrm{Fro}} ^ 2

and the regularized objective function is:

.. math::
    d_{\mathrm{Fro}}(X, WH)
    + \alpha \rho ||W||_1 + \alpha \rho ||H||_1
    + \frac{\alpha(1-\rho)}{2} ||W||_{\mathrm{Fro}} ^ 2
    + \frac{\alpha(1-\rho)}{2} ||H||_{\mathrm{Fro}} ^ 2

:class:`NMF` regularizes both W and H. The public function
:func:`non_negative_factorization` allows a finer control through the
:attr:`regularization` attribute, and may regularize only W, only H, or both.

NMF with a beta-divergence
--------------------------

As described previously, the most widely used distance function is the squared
Frobenius norm, which is an obvious extension of the Euclidean norm to
matrices:

.. math::
    d_{\mathrm{Fro}}(X, Y) = \frac{1}{2} ||X - Y||_{Fro}^2 = \frac{1}{2} \sum_{i,j} (X_{ij} - {Y}_{ij})^2

Other distance functions can be used in NMF as, for example, the (generalized)
Kullback-Leibler (KL) divergence, also referred as I-divergence:

.. math::
    d_{KL}(X, Y) = \sum_{i,j} (X_{ij} \log(\frac{X_{ij}}{Y_{ij}}) - X_{ij} + Y_{ij})

Or, the Itakura-Saito (IS) divergence:

.. math::
    d_{IS}(X, Y) = \sum_{i,j} (\frac{X_{ij}}{Y_{ij}} - \log(\frac{X_{ij}}{Y_{ij}}) - 1)

These three distances are special cases of the beta-divergence family, with
:math:`\beta = 2, 1, 0` respectively [6]_. The beta-divergence are
defined by :

.. math::
    d_{\beta}(X, Y) = \sum_{i,j} \frac{1}{\beta(\beta - 1)}(X_{ij}^\beta + (\beta-1)Y_{ij}^\beta - \beta X_{ij} Y_{ij}^{\beta - 1})

.. figure:: ../auto_examples/decomposition/images/sphx_glr_plot_beta_divergence_001.png
    :target: ../auto_examples/decomposition/plot_beta_divergence.html
    :align: center
    :scale: 75%

Note that this definition is not valid if :math:`\beta \in (0; 1)`, yet it can
be continuously extended to the definitions of :math:`d_{KL}` and :math:`d_{IS}`
respectively.

:class:`NMF` implements two solvers, using Coordinate Descent ('cd') [5]_, and
Multiplicative Update ('mu') [6]_. The 'mu' solver can optimize every
beta-divergence, including of course the Frobenius norm (:math:`\beta=2`), the
(generalized) Kullback-Leibler divergence (:math:`\beta=1`) and the
Itakura-Saito divergence (:math:`\beta=0`). Note that for
:math:`\beta \in (1; 2)`, the 'mu' solver is significantly faster than for other
values of :math:`\beta`. Note also that with a negative (or 0, i.e.
'itakura-saito') :math:`\beta`, the input matrix cannot contain zero values.

The 'cd' solver can only optimize the Frobenius norm. Due to the
underlying non-convexity of NMF, the different solvers may converge to
different minima, even when optimizing the same distance function.

NMF is best used with the ``fit_transform`` method, which returns the matrix W.
The matrix H is stored into the fitted model in the ``components_`` attribute;
the method ``transform`` will decompose a new matrix X_new based on these
stored components::

    >>> import numpy as np
    >>> X = np.array([[1, 1], [2, 1], [3, 1.2], [4, 1], [5, 0.8], [6, 1]])
    >>> from sklearn.decomposition import NMF
    >>> model = NMF(n_components=2, init='random', random_state=0)
    >>> W = model.fit_transform(X)
    >>> H = model.components_
    >>> X_new = np.array([[1, 0], [1, 6.1], [1, 0], [1, 4], [3.2, 1], [0, 4]])
    >>> W_new = model.transform(X_new)

.. topic:: Examples:

    * :ref:`sphx_glr_auto_examples_decomposition_plot_faces_decomposition.py`
    * :ref:`sphx_glr_auto_examples_applications_plot_topics_extraction_with_nmf_lda.py`
    * :ref:`sphx_glr_auto_examples_decomposition_plot_beta_divergence.py`

.. topic:: References:

    .. [1] `"Learning the parts of objects by non-negative matrix factorization"
      <http://www.columbia.edu/~jwp2128/Teaching/E4903/papers/nmf_nature.pdf>`_
      D. Lee, S. Seung, 1999

    .. [2] `"Non-negative Matrix Factorization with Sparseness Constraints"
      <http://www.jmlr.org/papers/volume5/hoyer04a/hoyer04a.pdf>`_
      P. Hoyer, 2004

    .. [4] `"SVD based initialization: A head start for nonnegative
      matrix factorization"
      <http://scgroup.hpclab.ceid.upatras.gr/faculty/stratis/Papers/HPCLAB020107.pdf>`_
      C. Boutsidis, E. Gallopoulos, 2008

    .. [5] `"Fast local algorithms for large scale nonnegative matrix and tensor
      factorizations."
      <http://www.bsp.brain.riken.jp/publications/2009/Cichocki-Phan-IEICE_col.pdf>`_
      A. Cichocki, A. Phan, 2009

    .. [6] `"Algorithms for nonnegative matrix factorization with the beta-divergence"
      <https://arxiv.org/pdf/1010.1763.pdf>`_
      C. Fevotte, J. Idier, 2011


.. _LatentDirichletAllocation:

Latent Dirichlet Allocation (LDA)
=================================

Latent Dirichlet Allocation is a generative probabilistic model for collections of
discrete dataset such as text corpora. It is also a topic model that is used for
discovering abstract topics from a collection of documents.

The graphical model of LDA is a three-level generative model:

.. image:: ../images/lda_model_graph.png
   :align: center

Note on notations presented in the graphical model above, which can be found in 
Hoffman et al. (2013):

  * The corpus is a collection of :math:`D` documents.
  * A document is a sequence of :math:`N` words.
  * There are :math:`K` topics in the corpus. 
  * The boxes represent repeated sampling. 

In the graphical model, each node is a random variable and has a role in the 
generative process. A shaded node indicates an observed variable and an unshaded 
node indicates a hidden (latent) variable. In this case, words in the corpus are 
the only data that we observe. The latent variables determine the random mixture 
of topics in the corpus and the distribution of words in the documents. 
The goal of LDA is to use the observed words to infer the hidden topic 
structure. 

When modeling text corpora, the model assumes the following generative process 
for a corpus with :math:`D` documents and :math:`K` topics, with :math:`K` 
corresponding to :attr:`n_components` in the API:

  1. For each topic :math:`k \in K`, draw :math:`\beta_k \sim 
     \mathrm{Dirichlet}(\eta)`. This provides a distribution over the words, 
     i.e. the probability of a word appearing in topic :math:`k`. 
     :math:`\eta` corresponds to :attr:`topic_word_prior`. 

  2. For each document :math:`d \in D`, draw the topic proportions 
     :math:`\theta_d \sim \mathrm{Dirichlet}(\alpha)`. :math:`\alpha` 
     corresponds to :attr:`doc_topic_prior`. 

  3. For each word :math:`i` in document :math:`d`:

    a. Draw the topic assignment :math:`z_{di} \sim \mathrm{Multinomial}
       (\theta_d)`
    b. Draw the observed word :math:`w_{ij} \sim \mathrm{Multinomial}
       (\beta_{z_{di}})`

For parameter estimation, the posterior distribution is:

.. math::
  p(z, \theta, \beta |w, \alpha, \eta) =
    \frac{p(z, \theta, \beta|\alpha, \eta)}{p(w|\alpha, \eta)}

Since the posterior is intractable, variational Bayesian method
uses a simpler distribution :math:`q(z,\theta,\beta | \lambda, \phi, \gamma)`
to approximate it, and those variational parameters :math:`\lambda`, 
:math:`\phi`, :math:`\gamma` are optimized to maximize the Evidence 
Lower Bound (ELBO):

.. math::
  \log\: P(w | \alpha, \eta) \geq L(w,\phi,\gamma,\lambda) \overset{\triangle}{=}
    E_{q}[\log\:p(w,z,\theta,\beta|\alpha,\eta)] - E_{q}[\log\:q(z, \theta, \beta)]

Maximizing ELBO is equivalent to minimizing the Kullback-Leibler(KL) divergence
between :math:`q(z,\theta,\beta)` and the true posterior
:math:`p(z, \theta, \beta |w, \alpha, \eta)`.

:class:`LatentDirichletAllocation` implements the online variational Bayes 
algorithm and supports both online and batch update methods.
While the batch method updates variational variables after each full pass through 
the data, the online method updates variational variables from mini-batch data 
points.

.. note::

  Although the online method is guaranteed to converge to a local optimum point, the quality of
  the optimum point and the speed of convergence may depend on mini-batch size and
  attributes related to learning rate setting.

When :class:`LatentDirichletAllocation` is applied on a "document-term" matrix, the matrix
will be decomposed into a "topic-term" matrix and a "document-topic" matrix. While
"topic-term" matrix is stored as :attr:`components_` in the model, "document-topic" matrix
can be calculated from ``transform`` method.

:class:`LatentDirichletAllocation` also implements ``partial_fit`` method. This is used
when data can be fetched sequentially.

.. topic:: Examples:

    * :ref:`sphx_glr_auto_examples_applications_plot_topics_extraction_with_nmf_lda.py`

.. topic:: References:

    * `"Latent Dirichlet Allocation"
      <http://www.jmlr.org/papers/volume3/blei03a/blei03a.pdf>`_
      D. Blei, A. Ng, M. Jordan, 2003

    * `"Online Learning for Latent Dirichlet Allocation”
      <https://papers.nips.cc/paper/3902-online-learning-for-latent-dirichlet-allocation.pdf>`_
      M. Hoffman, D. Blei, F. Bach, 2010

    * `"Stochastic Variational Inference"
      <http://www.columbia.edu/~jwp2128/Papers/HoffmanBleiWangPaisley2013.pdf>`_
      M. Hoffman, D. Blei, C. Wang, J. Paisley, 2013


See also :ref:`nca_dim_reduction` for dimensionality reduction with
Neighborhood Components Analysis.
