# %%
"""
============================
Image to Graph Conversion
============================

This example demonstrates how to convert a 2D image array into a graph
representation using the `img_to_graph` function from the
`sklearn.feature_extraction.image` module.
The graph is then normalized and the original image is displayed.

"""
# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause
# %%
import matplotlib.pyplot as plt
import numpy as np

from sklearn.feature_extraction import img_to_graph
from sklearn.preprocessing import normalize

# Create a sample image (2D array)
image = np.array([[0, 1, 2], [3, 4, 5], [6, 7, 8]])

# Convert the image to a graph
graph = img_to_graph(image)

# Normalize the graph
graph = normalize(graph, norm="l1", axis=1)

# Display the original image
plt.imshow(image, cmap="gray")
plt.title("Original Image")
plt.show()

# Print the graph representation
print("Graph representation (sparse matrix):")
print(graph)
