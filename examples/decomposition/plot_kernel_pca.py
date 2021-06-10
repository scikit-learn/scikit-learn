"""
==========
Kernel PCA
==========

This example shows the difference between the Principal Components Analysis
(:class:`~sklearn.decomposition.PCA`) and its kernalize version
(:class:`~sklearn.decomposition.KernelPCA`).

On the one hand, we show that :class:`~sklearn.decomposition.KernelPCA` is able
to find a projection of the data that makes them linearly separable while it is
not the case with :class:`~sklearn.decomposition.PCA`.

On the other hand, we show that inverting this projection is an
approximation with  :class:`~sklearn.decomposition.KernelPCA`,
while being exact with :class:`~sklearn.decomposition.PCA`.
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
# In this section, we will show the advantages of using a kernel when
# projecting data using a Principal Component Analysis (PCA). We create a
# dataset made of two nested circles.
from sklearn.datasets import make_circles
from sklearn.model_selection import train_test_split

X, y = make_circles(n_samples=1_000, factor=.3, noise=.05, random_state=0)
X_train, X_test, y_train, y_test = train_test_split(
    X, y, stratify=y, random_state=0)

# %%
# Let's have a quick first look to the dataset generated.
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
# The samples from each class cannot be linearly separated: we come with a
# straightline that would split the samples from the inner circle to outer
# circle. Perfectly, a potential decision function would be a circle separating
# both circles.
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

fig.subplots_adjust(wspace=0.3)

# %%
# We recall that PCA will project the data using a linear projection.
# Intuitively, it means that the coordinate system will be rotated with an some
# rescaling of the axis. This rescaling will depend on the variance of the
# data.
#
# Thus, looking at the projection made using PCA (i.e. the middle figure), we
# see that there is no change regarding the scaling; indeed the data being two
# concentric circles centered in zero, the variance of the original data was
# already maximized. However, we can see that the data have been rotated. As a
# conclusion, we see that such a projection would not help if define a linear
# classifier to distinguish samples from both classes.
#
# Using a kernel allows to make a non-linear projection. Here, by using an RBF
# kernel, we expect that the projection to unfold the dataset but keeping that
# point close in the original space should still be close in the new space.
#
# We observe such behaviour in the figure on the right: the samples of a given
# class are closer to each other than the samples from the opposite class. The
# "radial" effect make that we unrolled the circle. Now, we can use a linear
# classifier to separate the samples from the two classes.
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
# While we see a perfect reconstruction, we observe a different results for
# :class:`~sklearn.decomposition.KernelPCA`. Indeed,
# :meth:`~sklearn.decomposition.KernelPCA.inverse_transform` cannot rely on an
# analytical back-projection and thus an extact reconstruction. Instead, a
# :class:`~sklearn.linear_model.KernelRidge` was trained to learn a projection
# function to map a sample from the PCA basis into the original feature space.
# This method is therefore an approximation leading to small difference. The
# parameter `alpha` in the :class:`~sklearn.decomposition.KernelPCA` is used
# to penalized the mapping function to fit more or less the training data.
