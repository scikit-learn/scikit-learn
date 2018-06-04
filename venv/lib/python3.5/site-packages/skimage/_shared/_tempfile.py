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
    ...     im = np.zeros((5, 5), np.uint8)
    ...     io.imsave(tempfile, im)
    ...     assert np.all(io.imread(tempfile) == im)
    """
    tempfile_stream = NamedTemporaryFile(suffix=suffix, delete=False)
    tempfile = tempfile_stream.name
    tempfile_stream.close()
    yield tempfile
    os.remove(tempfile)
