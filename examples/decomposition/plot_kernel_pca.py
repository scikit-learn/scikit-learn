"""
==========
Kernel PCA
==========

This example shows the difference between Principal Components Analysis
(:class:`~sklearn.decomposition.PCA`) and
:class:`~sklearn.decomposition.KernelPCA`.

On the one hand, we show that :class:`~sklearn.decomposition.KernelPCA` is able
to find a projection of the data that makes them linearly separable while it is
not the case with
:class:`~sklearn.decomposition.PCA`.

On the other hand, we show that inverting this projection is an
approximation with  :class:`~sklearn.decomposition.KernelPCA`,
while being exact with :class:`~sklearn.decomposition.PCA`.

Finally, we show that this limitation can be useful in some applications such
as image denoising.
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
fig, axs = plt.subplots(ncols=2, nrows=2, figsize=(8, 8))

axs[0, 0].scatter(X_train[:, 0], X_train[:, 1], c=y_train)
axs[0, 0].set_ylabel("Feature #1")
axs[0, 0].set_xlabel("Feature #0")
axs[0, 0].set_title("Training data")

axs[0, 1].scatter(X_test[:, 0], X_test[:, 1], c=y_test)
axs[0, 1].set_xlabel("Feature #0")
axs[0, 1].set_title("Testing data")

axs[1, 0].scatter(X_test_pca[:, 0], X_test_pca[:, 1], c=y_test)
axs[1, 0].set_ylabel("Principal component #1")
axs[1, 0].set_xlabel("Principal component #0")
axs[1, 0].set_title("Projection of testing data\n using PCA")

axs[1, 1].scatter(X_test_kernel_pca[:, 0], X_test_kernel_pca[:, 1], c=y_test)
axs[1, 1].set_xlabel("Principal component #0")
axs[1, 1].set_title("Projection of testing data\n using KernelPCA")

fig.subplots_adjust(hspace=0.4)

# %%
# We recall that PCA will project the data using a linear projection.
# Intuitively, it means that the coordinate system will be rotated with an some
# rescaling of the axis. This rescaling will depend on the variance of the
# data.
#
# Thus, looking at the projection made using PCA (i.e. figure on the
# bottom-left), we see that there is no change regarding the scaling; indeed
# the data being two concentric circles centered in zero, the variance of the
# original data was already maximized. However, we can see that the data have
# been rotated. As a conclusion, we see that such a projection would not help
# define a linear classifier to distinguish samples from both classes.
#
# Using a kernel allows to make a non-linear projection. Here, by using an RBF
# kernel, we expect that the projection to unfold the dataset but keeping that
# point close in the original space should still be close in the new space.
#
# We observe such behaviour in the bottom-right figure: the samples of a given
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
#
# Application to image denoising
# ------------------------------
#
# In this section, we will show how one can use the approximation function
# learned to denoise image.

# %%
import numpy as np
from sklearn.datasets import fetch_openml
from sklearn.preprocessing import MinMaxScaler
from sklearn.model_selection import train_test_split

X, y = fetch_openml(data_id=41082, as_frame=False, return_X_y=True)
X_train, X_test, y_train, y_test = train_test_split(
    X, y, stratify=y, random_state=0, train_size=1_000, test_size=100
)
min_max_scaler = MinMaxScaler()
X_train = min_max_scaler.fit_transform(X_train)
X_test = min_max_scaler.transform(X_test)

rng = np.random.RandomState(0)
noise = rng.normal(scale=0.25, size=X_test.shape)
X_test_noisy = X_test + noise


# %%
def plot_digits(X, title):
    """Small helper function to plot 100 digits."""
    fig, axs = plt.subplots(nrows=10, ncols=10, figsize=(8, 8))
    for img, ax in zip(X, axs.ravel()):
        ax.imshow(img.reshape((16, 16)), cmap="Greys")
        ax.axis("off")
    fig.suptitle(title, fontsize=30)


# %%
plot_digits(X_train, "Uncorrupted train images")
plot_digits(X_test, "Uncorrupted test images")
plot_digits(X_test_noisy,
            f"Noisy test images - "
            f"MSE: {np.mean((X_test - X_test_noisy) ** 2):.2f}")

# %%
# We created a training and testing of 1,000 samples and a test set of 100
# samples. Also, we created a corrupted testing set that correspond to the
# original test set with additional Gaussian noise.
#
# The idea of this section, is to show that we can denoise the corrupted test
# set by a learning a PCA basis on the uncorrupted train set. We will use
# both a PCA and a kernel-based PCA.
pca = PCA(n_components=32)
kernel_pca = KernelPCA(n_components=200, kernel="rbf", gamma=1e-3,
                       fit_inverse_transform=True, alpha=5e-3)

pca.fit(X_train)
_ = kernel_pca.fit(X_train)

# %%
# Now, can transform and reconstruct the noisy test set. Since we used less
# components than the number of original features, we will get an approximation
# of the original set. Indeed, by dropping the components explaining less
# variance in PCA, we hope to remove noise. Similar thinking happen in kernel
# PCA; however, we expect a better reconstruction because we use a non-linear
# kernel to learn the PCA basis and a kernel ridge to learn the mapping
# function.
X_reconstructed_kernel_pca = kernel_pca.inverse_transform(
    kernel_pca.transform(X_test_noisy))
X_reconstructed_pca = pca.inverse_transform(pca.transform(X_test_noisy))

# %%
plot_digits(X_test, "Uncorrupted test images")
plot_digits(X_reconstructed_pca,
            f"PCA reconstruction - "
            f"MSE: {np.mean((X_test - X_reconstructed_pca) ** 2):.2f}")
plot_digits(X_reconstructed_kernel_pca,
            f"Kernel PCA reconstruction - "
            f"MSE: {np.mean((X_test - X_reconstructed_kernel_pca) ** 2):.2f}")

# %%
# Even if both PCA and kernel PCA have the same MSE, a qualitative analysis
# will favor the output of the kernel PCA. However, it should be noted that
# the results of the denoising with kernel PCA will depend of the parameters
# `n_components`, `gamma`, and `alpha`.
