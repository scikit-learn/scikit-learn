
from matplotlib.testing.decorators import image_comparison
from mpl_toolkits.axes_grid1 import ImageGrid
import numpy as np
import matplotlib.pyplot as plt


@image_comparison(baseline_images=['imagegrid_cbar_mode'],
                  extensions=['png'],
                  remove_text=True,
                  style='mpl20')
def test_imagegrid_cbar_mode_edge():
    X, Y = np.meshgrid(np.linspace(0, 6, 30), np.linspace(0, 6, 30))
    arr = np.sin(X) * np.cos(Y) + 1j*(np.sin(3*Y) * np.cos(Y/2.))

    fig = plt.figure(figsize=(18, 9))

    positions = (241, 242, 243, 244, 245, 246, 247, 248)
    directions = ['row']*4 + ['column']*4
    cbar_locations = ['left', 'right', 'top', 'bottom']*2

    for position, direction, location in zip(positions,
                                             directions,
                                             cbar_locations):
        grid = ImageGrid(fig, position,
                         nrows_ncols=(2, 2),
                         direction=direction,
                         cbar_location=location,
                         cbar_size='20%',
                         cbar_mode='edge')
        ax1, ax2, ax3, ax4, = grid

        im1 = ax1.imshow(arr.real, cmap='nipy_spectral')
        im2 = ax2.imshow(arr.imag, cmap='hot')
        im3 = ax3.imshow(np.abs(arr), cmap='jet')
        im4 = ax4.imshow(np.arctan2(arr.imag, arr.real), cmap='hsv')

        # Some of these colorbars will be overridden by later ones,
        # depending on the direction and cbar_location
        ax1.cax.colorbar(im1)
        ax2.cax.colorbar(im2)
        ax3.cax.colorbar(im3)
        ax4.cax.colorbar(im4)
