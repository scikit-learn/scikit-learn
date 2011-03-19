from .base import load_diabetes
from .base import load_digits
from .base import load_filenames
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


# backward compatibility
def load_files(*args, **kwargs):
    """Deprecated, use ``scikits.learn.datasets.load_filenames`` instead"""
    import warnings
    warnings.warn("Deprecated in 0.8 to be removed in 0.9: use load_filenames")
    return load_filenames(*args, **kwargs)
