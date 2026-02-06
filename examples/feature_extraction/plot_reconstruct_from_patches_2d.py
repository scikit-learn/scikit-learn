"""
============================
Reconstruct from Patches 2D Example
============================

This example demonstrates how to use the `reconstruct_from_patches_2d` function from the
`sklearn.feature_extraction.image` module to reconstruct an image from its patches.

"""

# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause
# %%
import matplotlib.pyplot as plt
import numpy as np

from sklearn.feature_extraction import image

# Create a sample image (2D array)
original_image = np.array([[0, 1, 2], [3, 4, 5], [6, 7, 8]])

# Extract patches from the image
patches = image.extract_patches_2d(original_image, (2, 2))

# Reconstruct the image from patches
reconstructed_image = image.reconstruct_from_patches_2d(patches, original_image.shape)

# Display the original and reconstructed images
fig, axes = plt.subplots(1, 2, figsize=(8, 4))
axes[0].imshow(original_image, cmap="gray")
axes[0].set_title("Original Image")
axes[1].imshow(reconstructed_image, cmap="gray")
axes[1].set_title("Reconstructed Image")
plt.show()

# Print the original and reconstructed images
print("Original Image:\n", original_image)
print("Reconstructed Image:\n", reconstructed_image)
