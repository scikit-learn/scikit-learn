"""
===========================
Vector Quantization Example
===========================

This example shows how one can use :class:`~sklearn.preprocessing.KBinsDiscretizer`
to perform vector quantization on a set of toy image, the raccoon face.
"""

# Authors: Gael Varoquaux
#          Jaques Grobler
# License: BSD 3 clause

# %%
# Original image
# --------------
#
# We start by loading the raccoon face image from SciPy. We will additionally check
# a couple of information regarding the image, such as the shape and data type used
# to store the image.
#
# Note that depending of the SciPy version, we have to adapt the import since the
# function returning the image is not located in the same module. Also, SciPy >= 1.10
# will required the package `pooch` to be installed.
try:  # Scipy >= 1.10
    from scipy.datasets import face
except ImportError:
    from scipy.misc import face

raccoon_face = face(gray=True)

print(f"The dimension of the image is {raccoon_face.shape}")
print(f"The data used to encode the image is of type {raccoon_face.dtype}")
print(f"The number of bytes taken in RAM is {raccoon_face.nbytes}")

# %%
# Thus the image is a 2D array of 768 pixels in height and 1024 pixels in width. Each
# value is a 8-bit unsigned integer, which means that the image is encoded using 8
# bits per pixel. The total memory usage of the image is 786 kilobytes (1 bytes equal
# 8 bits).
#
# Using 8-bit unsigned integer means that the image is encoded using 256 different
# shades of gray, at most. We can check the distribution of these values.
import matplotlib.pyplot as plt

fig, ax = plt.subplots(ncols=2, figsize=(12, 4))

ax[0].imshow(raccoon_face, cmap=plt.cm.gray)
ax[0].axis("off")
ax[0].set_title("Rendering of the image")
ax[1].hist(raccoon_face.ravel(), bins=256)
ax[1].set_xlabel("Pixel value")
ax[1].set_ylabel("Number of pixels")
ax[1].set_title("Distribution of the pixel values")
_ = fig.suptitle("Original image of a raccoon face")

# %%
# Compression via vector quantization
# -----------------------------------
# The idea behind compression via vector quantization is to use a smaller number of
# gray levels to represent the image. For instance, we can use 8 values instead of 256
# values. Therefore, it means that we could efficiently use 1 bit instead of 8 bits to
# encode a single pixel and therefore reduce the memory usage by a factor of 8. We will
# later discuss about this memory usage.
#
# Encoding strategy
# """""""""""""""""
#
# The compression can be performed by using a
# :class:`~sklearn.preprocessing.KBinsDiscretizer`. We need to choose of a strategy
# to define the 8 values. The simplest strategy is to define the 8 values equally
# spaced that correspond to `strategy="uniform"`. From the previous histogram, we know
# that this strategy is certainly not optimal.

from sklearn.preprocessing import KBinsDiscretizer

n_bins = 8
encoder = KBinsDiscretizer(
    n_bins=n_bins, encode="ordinal", strategy="uniform", random_state=0
)
compressed_raccoon_uniform = encoder.fit_transform(raccoon_face.reshape(-1, 1)).reshape(
    raccoon_face.shape
)

fig, ax = plt.subplots(ncols=2, figsize=(12, 4))
ax[0].imshow(compressed_raccoon_uniform, cmap=plt.cm.gray)
ax[0].axis("off")
ax[0].set_title("Rendering of the image")
ax[1].hist(compressed_raccoon_uniform.ravel(), bins=256)
ax[1].set_xlabel("Pixel value")
ax[1].set_ylabel("Number of pixels")
ax[1].set_title("Distribution of the pixel values")
_ = fig.suptitle("Raccoon face compressed using 1-bit and a uniform strategy")

# %%
# Qualitatively, we can spot some small regions where we see the effect of the
# compression (e.g. leaves on the bottom right corner). But after all, the resulting
# image is still looking good.
#
# Looking a the distribution, we observe that all pixels values have been mapped to
# 8 different values. But we can also see the amount of pixel in each bin is different.
# This might not be the most effective strategy. We can check the correspondance between
# the 8 values and the original pixel values.
bin_edges = encoder.bin_edges_[0]
bin_center = bin_edges[:-1] + (bin_edges[1:] - bin_edges[:-1]) / 2
bin_center

# %%
plt.hist(raccoon_face.ravel(), bins=256)
for center in bin_center:
    plt.axvline(center, color="red")

# %%
# As previously stated, the strategy used is not optimal. The value 7, will encode a
# rather small amount of pixels while the value 3 will encode a large amount of pixels.
# We can instead use a clustering strategy such as k-means to find the 8 values.
encoder = KBinsDiscretizer(
    n_bins=n_bins, encode="ordinal", strategy="kmeans", random_state=0
)
compressed_raccoon_kmeans = encoder.fit_transform(raccoon_face.reshape(-1, 1)).reshape(
    raccoon_face.shape
)

fig, ax = plt.subplots(ncols=2, figsize=(12, 4))
ax[0].imshow(compressed_raccoon_kmeans, cmap=plt.cm.gray)
ax[0].axis("off")
ax[0].set_title("Rendering of the image")
ax[1].hist(compressed_raccoon_kmeans.ravel(), bins=256)
ax[1].set_xlabel("Pixel value")
ax[1].set_ylabel("Number of pixels")
ax[1].set_title("Distribution of the pixel values")
_ = fig.suptitle("Raccoon face compressed using 1-bit and a K-means strategy")

# %%
bin_edges = encoder.bin_edges_[0]
bin_center = bin_edges[:-1] + (bin_edges[1:] - bin_edges[:-1]) / 2
bin_center

# %%
plt.hist(raccoon_face.ravel(), bins=256)
for center in bin_center:
    plt.axvline(center, color="red")

# %%
# The bins are now more balanced and the center of each bins are different from
# the uniform strategy. Note that we could enforce the same number of pixels
# per bin by using the `strategy="quantile"` instead of `strategy="kmeans"`.
#
# Memory footprint
# """"""""""""""""
#
# We can quickly discussed the memory usage of the compressed image. We previously
# stated that we should save 8 times less memory. Let's check it.
print(f"The number of bytes taken in RAM is {compressed_raccoon_kmeans.nbytes}")
print(f"Compression ratio: {compressed_raccoon_kmeans.nbytes / raccoon_face.nbytes}")

# %%
# This is quite surprising to see that our compressed image is taking x8 more memory
# than the original image. This is indeed the opposit of what we expected. The reason
# is mainly due to the type of data used to encode the image.
print(f"Type of the compressed image: {compressed_raccoon_kmeans.dtype}")

# %%
# Indeed, the output of the :class:`~sklearn.preprocessing.KBinsDiscretizer` is an
# array of 64-bit float. It means that it takes x8 more memory. However, we use this
# 64-bit float representation to encode 8 values. Indeed, we will save memory only if
# we cast the compressed image into an array of 1-bit integer. We could use potentially
# the function `as_type` from `numpy`. However, a 1-bit integer representation does not
# exist and to encode the 8 values, we will need to use the 8-bit unsigned integer
# representation as well.
#
# In reality, the above compression would make sense with 64-bit float representation
# of original image to observe a memory gain.
