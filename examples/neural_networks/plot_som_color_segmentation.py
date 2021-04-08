"""
========================================================
Example of color quantization of an image using SOM
========================================================

In this example, we use a rectangular SOM with 6 x 6 to 
learn all mojority colors of an image. The example try to
classify each pixel of the image into each neuron mapped
into this SOM network. The result is a segmented image with
all 36 predominant colors.
"""
print(__doc__)

# Authors: Julio Faracco
# License: BSD

import numpy as np
import matplotlib.pyplot as plt

from skimage import io

from sklearn.neural_network import SOM

image_url = "https://upload.wikimedia.org/wikipedia/commons/thumb/8/8f/Bachalpsee_reflection.jpg/250px-Bachalpsee_reflection.jpg"

image_np = io.imread(image_url)
x, y, input_len = image_np.shape

# reshaping the pixels matrix
pixels = np.reshape(image_np, (x * y, input_len)) / 255.

print('Training SOM...')
som = SOM(6, 6, 3, sigma=1., verbose=0,
          learning_rate=0.2, neighborhood_function='bubble', random_state=42)

som.randomize_weights_from_data(pixels)
starting_weights = som.weights().copy()
som.fit(pixels)

print('Building new image...')
clustered = som.predict(pixels)
clustered = np.reshape(clustered, image_np.shape)
print('Done.')

# show the result
plt.figure(figsize=(7, 7))
plt.figure(1)
plt.subplot(221)
plt.title('Original')
plt.imshow(image_np)
plt.subplot(222)
plt.title('Result')
plt.imshow(clustered)

plt.subplot(223)
plt.title('Initial colors')
plt.imshow(starting_weights, interpolation='none')
plt.subplot(224)
plt.title('learned colors')
plt.imshow(som.weights(), interpolation='none')

plt.tight_layout()
plt.show()
