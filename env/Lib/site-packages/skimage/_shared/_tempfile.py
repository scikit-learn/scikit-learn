from tempfile import NamedTemporaryFile
from contextlib import contextmanager
import os


@contextmanager
def temporary_file(suffix=''):
    """Yield a writeable temporary filename that is deleted on context exit.

    Parameters
    ----------
    suffix : string, optional
        The suffix for the file.

    Examples
    --------
    >>> import numpy as np
    >>> from skimage import io
    >>> with temporary_file('.tif') as tempfile:
    ...     im = np.arange(25, dtype=np.uint8).reshape((5, 5))
    ...     io.imsave(tempfile, im)
    ...     assert np.all(io.imread(tempfile) == im)
    """
    with NamedTemporaryFile(suffix=suffix, delete=False) as tempfile_stream:
        tempfile = tempfile_stream.name

    yield tempfile
    os.remove(tempfile)
