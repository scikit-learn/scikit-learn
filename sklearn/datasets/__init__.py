"""
The :mod:`sklearn.datasets` module includes utilities to load datasets,
including methods to load and fetch popular reference datasets. It also
features some artificial data generators.
"""

from .base import load_diabetes
from .base import load_digits
from .base import load_files
from .base import load_iris
from .base import load_linnerud
from .base import load_boston
from .base import get_data_home
from .base import clear_data_home
from .base import load_sample_images
from .base import load_sample_image
from .mlcomp import load_mlcomp
from .lfw import load_lfw_pairs
from .lfw import load_lfw_people
from .lfw import fetch_lfw_pairs
from .lfw import fetch_lfw_people
from .twenty_newsgroups import fetch_20newsgroups
from .twenty_newsgroups import fetch_20newsgroups_vectorized
from .twenty_newsgroups import load_20newsgroups
from .mldata import fetch_mldata, mldata_filename
from .samples_generator import make_classification
from .samples_generator import make_multilabel_classification
from .samples_generator import make_regression
from .samples_generator import make_blobs
from .samples_generator import make_friedman1
from .samples_generator import make_friedman2
from .samples_generator import make_friedman3
from .samples_generator import make_low_rank_matrix
from .samples_generator import make_sparse_coded_signal
from .samples_generator import make_sparse_uncorrelated
from .samples_generator import make_spd_matrix
from .samples_generator import make_swiss_roll
from .samples_generator import make_s_curve
from .samples_generator import make_sparse_spd_matrix
from .svmlight_format import load_svmlight_file
from .svmlight_format import load_svmlight_files
from .svmlight_format import dump_svmlight_file
from .olivetti_faces import fetch_olivetti_faces
from .species_distributions import fetch_species_distributions
from ..utils import deprecated


# backward compatibility
@deprecated("to be removed in 0.9;"
            " use sklearn.datasets.load_files instead")
def load_filenames(*args, **kwargs):
    """Deprecated, use ``sklearn.datasets.load_files`` instead"""
    return load_files(*args, **kwargs)
