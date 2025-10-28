"""Image Processing for Python

scikit-image (a.k.a. ``skimage``) is a collection of algorithms for image
processing and computer vision.

Attributes
----------
__version__ : str
    The scikit-image version string.

Subpackages
-----------
color
    Color space conversion.
data
    Example images and datasets.
draw
    Drawing primitives, such as lines, circles, text, etc.
exposure
    Image intensity adjustment, e.g., histogram equalization, etc.
feature
    Feature detection and extraction, e.g., texture analysis, corners, etc.
filters
    Sharpening, edge finding, rank filters, thresholding, etc.
future
    Functionality with an experimental API.
graph
    Graph-based operations, e.g., shortest paths.
io
    Reading and saving of images and videos.
measure
    Measurement of image properties, e.g., region properties, contours.
metrics
    Metrics corresponding to images, e.g., distance metrics, similarity, etc.
morphology
    Morphological algorithms, e.g., closing, opening, skeletonization.
registration
    Image registration algorithms, e.g., optical flow or phase cross correlation.
restoration
    Restoration algorithms, e.g., deconvolution algorithms, denoising, etc.
segmentation
    Algorithms to partition images into meaningful regions or boundaries.
transform
    Geometric and other transformations, e.g., rotations, Radon transform.
util
    Generic utilities.
"""

__version__ = '0.25.2'

import lazy_loader as _lazy

__getattr__, *_ = _lazy.attach_stub(__name__, __file__)


# Don't use the `__all__` and `__dir__` returned by `attach_stubs` since that
# one would expose utility functions we don't want to advertise in our
# top-level module anymore.
__all__ = [
    "__version__",
    "color",
    "data",
    "draw",
    "exposure",
    "feature",
    "filters",
    "future",
    "graph",
    "io",
    "measure",
    "metrics",
    "morphology",
    "registration",
    "restoration",
    "segmentation",
    "transform",
    "util",
]


def __dir__():
    return __all__.copy()


# Logic for checking for improper install and importing while in the source
# tree when package has not been installed inplace.
# Code adapted from scikit-learn's __check_build module.
_INPLACE_MSG = """
It appears that you are importing a local scikit-image source tree. For
this, you need to have an inplace install. Maybe you are in the source
directory and you need to try from another location."""

_STANDARD_MSG = """
Your install of scikit-image appears to be broken.
Try re-installing the package following the instructions at:
https://scikit-image.org/docs/stable/user_guide/install.html"""


def _raise_build_error(e):
    # Raise a comprehensible error
    import os.path as osp

    local_dir = osp.split(__file__)[0]
    msg = _STANDARD_MSG
    if local_dir == "skimage":
        # Picking up the local install: this will work only if the
        # install is an 'inplace build'
        msg = _INPLACE_MSG
    raise ImportError(
        f"{e}\nIt seems that scikit-image has not been built correctly.\n{msg}"
    )


def _try_append_commit_info(version):
    """Append last commit date and hash to `version`, if available."""
    import subprocess
    from pathlib import Path

    try:
        output = subprocess.check_output(
            ['git', 'log', '-1', '--format="%h %aI"'],
            cwd=Path(__file__).parent,
            text=True,
        )
        if output:
            git_hash, git_date = (
                output.strip().replace('"', '').split('T')[0].replace('-', '').split()
            )
            version = '+'.join(
                [tag for tag in version.split('+') if not tag.startswith('git')]
            )
            version += f'+git{git_date}.{git_hash}'

    except (FileNotFoundError, subprocess.CalledProcessError):
        pass
    except OSError:
        pass  # If skimage is built with emscripten which does not support processes

    return version


if 'dev' in __version__:
    __version__ = _try_append_commit_info(__version__)


from skimage._shared.tester import PytestTester as _PytestTester

test = _PytestTester(__name__)
