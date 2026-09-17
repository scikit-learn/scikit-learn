"""
==============================================
Using multi-layer perceptron as an autoencoder
==============================================

This example shows how to use an :class:`~sklearn.neural_network.MLPRegressor`
as an autoencoder to denoise images. We will use the digits dataset for this
purpose.
"""

# %%
# The image denoising problem
# ---------------------------
#
# We use the digits dataset to illustrate the image denoising. We load the
# dataset using the :func:`sklearn.datasets.load_digits` function. We
# additionally scale the pixel intensities to the range [0, 1].
from sklearn.datasets import load_digits
from sklearn.preprocessing import minmax_scale

X, y = load_digits(return_X_y=True)
X = minmax_scale(X)

# %%
# We create a small helper function to plot the 40 first digits of a provided
# dataset.
import matplotlib.pyplot as plt


def plot_digits(X, y, title=None):
    fig, axes = plt.subplots(
        4, 10, figsize=(10, 4), subplot_kw={"xticks": (), "yticks": ()}
    )
    for i, ax in enumerate(axes.flat):
        ax.imshow(X[i].reshape(8, 8), cmap="binary", interpolation="nearest")
        ax.text(0.05, 0.05, str(y[i]), transform=ax.transAxes, color="green")

    if title:
        fig.suptitle(title)


# %%
# We plot the original digits dataset. We can see that the images are not
# corrupted.
plot_digits(X, y, title="Original digits dataset")

# %%
# To simulate a corrupted dataset, we add some noise to the original images.
import numpy as np

rng = np.random.RandomState(0)
X_noisy = X + rng.normal(scale=0.2, size=X.shape)
plot_digits(X_noisy, y, title="Noisy digits dataset")

# %%
# Before to design our autoencoder, we will first split the dataset into a
# training and a testing set to later evaluate the performance of this model.

from sklearn.model_selection import train_test_split

X_train, X_test, X_noisy_train, X_noisy_test, y_train, y_test = train_test_split(
    X, X_noisy, y, random_state=0
)

# %%
# Design the autoencoder
# ----------------------
#
# We define an autoencoder using a
# :class:`~sklearn.neural_network.MLPRegressor`. An autoencoder is a neural
# network that tries to reconstruct its input. In a a denoising application,
# the input is the noisy dataset and the output is the original uncorrupted
# dataset.
#
# Larger the number of layers and neurons, the more complex the model can be,
# and thus the best it can reconstruct the output. However, it will also be
# more prone to overfitting. Here, we will use 3 hidden layers with 300, 100
# and 300 neurons.
#
# We will `fit` the neural network to find the best weights to transform the
# noisy dataset into the original dataset. To evaluate our model, we will
# reconstruct the noisy test dataset.
from sklearn.neural_network import MLPRegressor

autoencoder = MLPRegressor(
    hidden_layer_sizes=[300, 100, 300],
    activation="relu",
    solver="adam",
    max_iter=1_000,
    random_state=0,
)

autoencoder.fit(X_noisy_train, X_train)
X_denoised_test = autoencoder.predict(X_noisy_test)

# %%
# To evaluate the performance of the model we will compute the peak
# signal-to-noise ratio (PSNR) which is a metric in decibels (dB). The higher
# the PSNR, the better the denoising.


def psnr_score(X_true, X_reconstructed):
    mse = np.mean((X_true - X_reconstructed) ** 2, axis=1)
    return 10 * np.log10(1 / mse)


# %%
# Now, we use the PSNR and plot the reconstructed image to evaluate the
# performance of the model quantitatively and qualitatively, respectively.
plot_digits(
    X_noisy_test,
    y_test,
    title=f"PSNR noisy dataset: {psnr_score(X_test, X_noisy_test).mean():.2f} dB",
)
plot_digits(
    X_denoised_test,
    y_test,
    title=f"PSNR denoised dataset: {psnr_score(X_test, X_denoised_test).mean():.2f} dB",
)

# %%
# We observe that the PSNR is higher for the denoised dataset and we confirm
# visually that the denoised dataset shows smoother images.
