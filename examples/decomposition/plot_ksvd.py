import matplotlib.pyplot as plt
import numpy as np
import scipy.misc

from sklearn.feature_extraction.image import (extract_patches_2d,
                                              reconstruct_from_patches_2d)
from sklearn.decomposition.ksvd import ksvd
from sklearn.linear_model import orthogonal_mp

plt.rc('image', interpolation='none', cmap='gray')


source_image = scipy.misc.lena() / 256.
plt.figure(figsize=(5, 5))
plt.imshow(source_image)


PATCH_SHAPE = (6, 6)


def show_atoms(atoms, rows, cols):
    assert rows * cols == atoms.shape[0]
    
    plt.figure(figsize=(cols, rows))
    for index in range(atoms.shape[0]):
        plt.subplot(rows, cols, index + 1)
        atom_image = atoms[index, :].reshape(PATCH_SHAPE)
        plt.imshow(atom_image)
        plt.xticks(())
        plt.yticks(())


patches_examples = extract_patches_2d(source_image, PATCH_SHAPE)[512*8:]
patches_examples = patches_examples.reshape(patches_examples.shape[0], -1)
print("Examples count:", patches_examples.shape[0])

train_sample = np.random.random_integers(0, patches_examples.shape[0], size=10000)
codes, dictionary = ksvd(patches_examples[train_sample, :], n_nonzero_coefs=8,
                         n_atoms=100, iteration_count=10, random_state=4132, n_jobs=4)
show_atoms(dictionary, 10, 10)


noisy_image = (source_image + 0.1 * np.random.randn(*source_image.shape))[256:-128, 256:-128]
plt.figure(figsize=(5, 5))
plt.imshow(noisy_image)


def restore(patch):
    patch_vector = patch.reshape(patch.size)
    code = orthogonal_mp(dictionary.T, patch_vector, n_nonzero_coefs=4)
    restored_vector = np.dot(dictionary.T, code)
    return restored_vector.reshape((1,) + PATCH_SHAPE)

def restore_image(image):
    patches = extract_patches_2d(image, PATCH_SHAPE)
    restored_patches = np.vstack([restore(patches[i, :, :]) for i in range(patches.shape[0])])
    return reconstruct_from_patches_2d(restored_patches, image.shape)

restored_image = restore_image(noisy_image)
plt.figure(figsize=(5, 5))
plt.imshow(restored_image)

