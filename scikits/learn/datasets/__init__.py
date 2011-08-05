from .base import load_diabetes
from .base import load_digits
from .base import load_files
from .base import load_iris
from .base import load_linnerud
from .base import get_data_home
from .base import clear_data_home
from .mlcomp import load_mlcomp
from .lfw import load_lfw_pairs
from .lfw import load_lfw_people
from .lfw import fetch_lfw_pairs
from .lfw import fetch_lfw_people
from .twenty_newsgroups import fetch_20newsgroups
from .twenty_newsgroups import load_20newsgroups
from .mldata import fetch_mldata, mldata_filename
from .samples_generator import swiss_roll
from .svmlight_format import load_svmlight_file
from .olivetti_faces import fetch_olivetti_faces
from ..utils import deprecated

# backward compatibility
@deprecated("to be removed in 0.9;"
            " use scikits.learn.datasets.load_files instead")
def load_filenames(*args, **kwargs):
    """Deprecated, use ``scikits.learn.datasets.load_files`` instead"""
    return load_files(*args, **kwargs)
