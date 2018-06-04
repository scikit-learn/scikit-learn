from __future__ import absolute_import, division, print_function

import six
from six import BytesIO
import glob
import os
import numpy as np
import pytest

from matplotlib.testing.decorators import image_comparison
from matplotlib import pyplot as plt
import matplotlib.cm as cm
import sys
on_win = (sys.platform == 'win32')


@image_comparison(baseline_images=['pngsuite'], extensions=['png'],
                  tol=0.03)
def test_pngsuite():
    dirname = os.path.join(
        os.path.dirname(__file__),
        'baseline_images',
        'pngsuite')
    files = sorted(glob.iglob(os.path.join(dirname, 'basn*.png')))

    fig = plt.figure(figsize=(len(files), 2))

    for i, fname in enumerate(files):
        data = plt.imread(fname)
        cmap = None  # use default colormap
        if data.ndim == 2:
            # keep grayscale images gray
            cmap = cm.gray
        plt.imshow(data, extent=[i, i + 1, 0, 1], cmap=cmap)

    plt.gca().patch.set_facecolor("#ddffff")
    plt.gca().set_xlim(0, len(files))


def test_imread_png_uint16():
    from matplotlib import _png
    img = _png.read_png_int(os.path.join(os.path.dirname(__file__),
                                     'baseline_images/test_png/uint16.png'))

    assert (img.dtype == np.uint16)
    assert np.sum(img.flatten()) == 134184960


def test_truncated_file(tmpdir):
    d = tmpdir.mkdir('test')
    fname = str(d.join('test.png'))
    fname_t = str(d.join('test_truncated.png'))
    plt.savefig(fname)
    with open(fname, 'rb') as fin:
        buf = fin.read()
    with open(fname_t, 'wb') as fout:
        fout.write(buf[:20])

    with pytest.raises(Exception):
        plt.imread(fname_t)


def test_truncated_buffer():
    b = BytesIO()
    plt.savefig(b)
    b.seek(0)
    b2 = BytesIO(b.read(20))
    b2.seek(0)

    with pytest.raises(Exception):
        plt.imread(b2)
