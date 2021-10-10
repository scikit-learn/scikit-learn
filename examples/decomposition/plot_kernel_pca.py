"""
==========
Kernel PCA
==========

This example shows the difference between the Principal Components Analysis
(:class:`~sklearn.decomposition.PCA`) and its kernalized version
(:class:`~sklearn.decomposition.KernelPCA`).

On the one hand, we show that :class:`~sklearn.decomposition.KernelPCA` is able
to find a projection of the data that makes them linearly separable while it is
not the case with :class:`~sklearn.decomposition.PCA`.

Finally, we show that inverting this projection is an approximation with
:class:`~sklearn.decomposition.KernelPCA`, while inverting is exact with
:class:`~sklearn.decomposition.PCA`.
"""
print(__doc__)

# Authors: Mathieu Blondel
#          Andreas Mueller
#          Guillaume Lemaitre
# License: BSD 3 clause

# %%
# Projecting data: `PCA` vs. `KernelPCA`
# --------------------------------------
#
# In this section, we show the advantages of using a kernel when
# projecting data using a Principal Component Analysis (PCA). We create a
# dataset made of two nested circles.
from sklearn.datasets import make_circles
from sklearn.model_selection import train_test_split

X, y = make_circles(n_samples=1_000, factor=.3, noise=.05, random_state=0)
X_train, X_test, y_train, y_test = train_test_split(
    X, y, stratify=y, random_state=0)

# %%
# Let's have a quick first look at the generated dataset.
import matplotlib.pyplot as plt

_, axs = plt.subplots(ncols=2, sharex=True, sharey=True, figsize=(8, 4))

axs[0].scatter(X_train[:, 0], X_train[:, 1], c=y_train)
axs[0].set_ylabel("Feature #1")
axs[0].set_xlabel("Feature #0")
axs[0].set_title("Training data")

axs[1].scatter(X_test[:, 0], X_test[:, 1], c=y_test)
axs[1].set_xlabel("Feature #0")
_ = axs[1].set_title("Testing data")

# %%
# The samples from each class cannot be linearly separated: there is no
# straight line that can split the samples of the inner set from the outer
# set.
#
# Now, we will use PCA with and without a kernel to see what is the effect of
# using such a kernel. The kernel used here is a radial basis function (RBF)
# kernel.
from sklearn.decomposition import PCA, KernelPCA

pca = PCA(n_components=2)
kernel_pca = KernelPCA(
    n_components=None, kernel="rbf", gamma=10, fit_inverse_transform=True,
    alpha=0.1)

X_test_pca = pca.fit(X_train).transform(X_test)
X_test_kernel_pca = kernel_pca.fit(X_train).transform(X_test)

# %%
fig, axs = plt.subplots(ncols=3, figsize=(14, 4))

axs[0].scatter(X_test[:, 0], X_test[:, 1], c=y_test)
axs[0].set_ylabel("Feature #1")
axs[0].set_xlabel("Feature #0")
axs[0].set_title("Testing data")

axs[1].scatter(X_test_pca[:, 0], X_test_pca[:, 1], c=y_test)
axs[1].set_ylabel("Principal component #1")
axs[1].set_xlabel("Principal component #0")
axs[1].set_title("Projection of testing data\n using PCA")

axs[2].scatter(X_test_kernel_pca[:, 0], X_test_kernel_pca[:, 1], c=y_test)
axs[2].set_ylabel("Principal component #1")
axs[2].set_xlabel("Principal component #0")
axs[2].set_title("Projection of testing data\n using KernelPCA")

# %%
# We recall that PCA projects the data linearly. Intuitively, it means that
# the coordinate system will be rotated after centering and rescaling on each
# axis. This rescaling will depend on the variance of the data.
#
# Thus, looking at the projection made using PCA (i.e. the middle figure), we
# see that there is no change regarding the scaling; indeed the data being two
# concentric circles centered in zero, the variance of the original data was
# already maximized. However, we can see that the data have been rotated. As a
# conclusion, we see that such a projection would not help if define a linear
# classifier to distinguish samples from both classes.
#
# Using a kernel allows to make a non-linear projection. Here, by using an RBF
# kernel, we expect that the projection will unfold the dataset while keeping
# approximately preserving the relative distances of pairs of data points that
# are close to one another in the original space.
#
# We observe such behaviour in the figure on the right: the samples of a given
# class are closer to each other than the samples from the opposite class,
# untangling both sample sets. Now, we can use a linear classifier to separate
# the samples from the two classes.
#
# Projecting into the original feature space
# ------------------------------------------
#
# One particularity to have in mind when using
# :class:`~sklearn.decomposition.KernelPCA` is related to the reconstruction
# (i.e. the back projection in the original feature space). With
# :class:`~sklearn.decomposition.PCA`, the reconstruction will be exact if
# `n_components` is the same than the number of original features as in this
# example. Thus, projecting the data on the PCA basis and projecting back will
# give the same dataset.
#
# We can investigate if we get a similar outcome with
# :class:`~sklearn.decomposition.KernelPCA`.
X_reconstructed_pca = pca.inverse_transform(pca.transform(X_test))
X_reconstructed_kernel_pca = kernel_pca.inverse_transform(
    kernel_pca.transform(X_test))

# %%
fig, axs = plt.subplots(ncols=3, sharex=True, sharey=True, figsize=(13, 4))

axs[0].scatter(X_test[:, 0], X_test[:, 1], c=y_test)
axs[0].set_ylabel("Feature #1")
axs[0].set_xlabel("Feature #0")
axs[0].set_title("Original test data")

axs[1].scatter(X_reconstructed_pca[:, 0], X_reconstructed_pca[:, 1], c=y_test)
axs[1].set_xlabel("Feature #0")
axs[1].set_title("Reconstruction via PCA")

axs[2].scatter(X_reconstructed_kernel_pca[:, 0],
               X_reconstructed_kernel_pca[:, 1], c=y_test)
axs[2].set_xlabel("Feature #0")
_ = axs[2].set_title("Reconstruction via KernelPCA")

# %%
# While we see a perfect reconstruction with
# :class:`~sklearn.decomposition.PCA` we observe a different results for
# :class:`~sklearn.decomposition.KernelPCA`.
#
# Indeed, :meth:`~sklearn.decomposition.KernelPCA.inverse_transform` cannot
# rely on an analytical back-projection and thus an extact reconstruction.
# Instead, a :class:`~sklearn.kernel_ridge.KernelRidge` is internally trained
# to learn a mapping from the PCA basis to the original feature space. This
# method is therefore an approximation introducing small differences when
# attempting to reconstruct the original input.
#
# To improve the reconstruction using
# :meth:`~sklearn.decomposition.KernelPCA.inverse_transform`, one can tune
# `alpha` in :class:`~sklearn.decomposition.KernelPCA`, the regularization term
# which controls the reliance on the training data during the training of
# the mapping.
