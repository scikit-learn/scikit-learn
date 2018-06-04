from matplotlib.backend_bases import FigureCanvasBase
from matplotlib.backend_bases import RendererBase

import matplotlib.pyplot as plt
import matplotlib.transforms as transforms
import matplotlib.path as path

import numpy as np
import os
import shutil
import tempfile


def test_uses_per_path():
    id = transforms.Affine2D()
    paths = [path.Path.unit_regular_polygon(i) for i in range(3, 7)]
    tforms = [id.rotate(i) for i in range(1, 5)]
    offsets = np.arange(20).reshape((10, 2))
    facecolors = ['red', 'green']
    edgecolors = ['red', 'green']

    def check(master_transform, paths, all_transforms,
              offsets, facecolors, edgecolors):
        rb = RendererBase()
        raw_paths = list(rb._iter_collection_raw_paths(
            master_transform, paths, all_transforms))
        gc = rb.new_gc()
        ids = [path_id for xo, yo, path_id, gc0, rgbFace in
               rb._iter_collection(gc, master_transform, all_transforms,
                                   range(len(raw_paths)), offsets,
                                   transforms.IdentityTransform(),
                                   facecolors, edgecolors, [], [], [False],
                                   [], 'data')]
        uses = rb._iter_collection_uses_per_path(
            paths, all_transforms, offsets, facecolors, edgecolors)
        if raw_paths:
            seen = np.bincount(ids, minlength=len(raw_paths))
            assert set(seen).issubset([uses - 1, uses])

    check(id, paths, tforms, offsets, facecolors, edgecolors)
    check(id, paths[0:1], tforms, offsets, facecolors, edgecolors)
    check(id, [], tforms, offsets, facecolors, edgecolors)
    check(id, paths, tforms[0:1], offsets, facecolors, edgecolors)
    check(id, paths, [], offsets, facecolors, edgecolors)
    for n in range(0, offsets.shape[0]):
        check(id, paths, tforms, offsets[0:n, :], facecolors, edgecolors)
    check(id, paths, tforms, offsets, [], edgecolors)
    check(id, paths, tforms, offsets, facecolors, [])
    check(id, paths, tforms, offsets, [], [])
    check(id, paths, tforms, offsets, facecolors[0:1], edgecolors)


def test_get_default_filename():
    try:
        test_dir = tempfile.mkdtemp()
        plt.rcParams['savefig.directory'] = test_dir
        fig = plt.figure()
        canvas = FigureCanvasBase(fig)
        filename = canvas.get_default_filename()
        assert filename == 'image.png'
    finally:
        shutil.rmtree(test_dir)


def test_get_default_filename_already_exists():
    # From #3068: Suggest non-existing default filename
    try:
        test_dir = tempfile.mkdtemp()
        plt.rcParams['savefig.directory'] = test_dir
        fig = plt.figure()
        canvas = FigureCanvasBase(fig)

        # create 'image.png' in figure's save dir
        open(os.path.join(test_dir, 'image.png'), 'w').close()

        filename = canvas.get_default_filename()
        assert filename == 'image-1.png'
    finally:
        shutil.rmtree(test_dir)
