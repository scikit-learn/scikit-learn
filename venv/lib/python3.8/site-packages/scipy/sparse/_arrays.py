from ._bsr import bsr_matrix
from ._coo import coo_matrix
from ._csc import csc_matrix
from ._csr import csr_matrix
from ._dia import dia_matrix
from ._dok import dok_matrix
from ._lil import lil_matrix


class _sparray:
    """This class provides a base class for all sparse arrays.

    It cannot be instantiated.  Most of the work is provided by subclasses.
    """

    _is_array = True

    @property
    def _bsr_container(self):
        return bsr_array

    @property
    def _coo_container(self):
        return coo_array

    @property
    def _csc_container(self):
        return csc_array

    @property
    def _csr_container(self):
        return csr_array

    @property
    def _dia_container(self):
        return dia_array

    @property
    def _dok_container(self):
        return dok_array

    @property
    def _lil_container(self):
        return lil_array

    # Restore elementwise multiplication
    def __mul__(self, *args, **kwargs):
        return self.multiply(*args, **kwargs)

    def __rmul__(self, *args, **kwargs):
        return self.multiply(*args, **kwargs)


def _matrix_doc_to_array(docstr):
    # For opimized builds with stripped docstrings
    if docstr is None:
        return None
    return docstr.replace("matrix", "array").replace("matrices", "arrays")


class bsr_array(_sparray, bsr_matrix):
    pass


class coo_array(_sparray, coo_matrix):
    pass


class csc_array(_sparray, csc_matrix):
    pass


class csr_array(_sparray, csr_matrix):
    pass


class dia_array(_sparray, dia_matrix):
    pass


class dok_array(_sparray, dok_matrix):
    pass


class lil_array(_sparray, lil_matrix):
    pass


bsr_array.__doc__ = _matrix_doc_to_array(bsr_matrix.__doc__)
coo_array.__doc__ = _matrix_doc_to_array(coo_matrix.__doc__)
csc_array.__doc__ = _matrix_doc_to_array(csc_matrix.__doc__)
csr_array.__doc__ = _matrix_doc_to_array(csr_matrix.__doc__)
dia_array.__doc__ = _matrix_doc_to_array(dia_matrix.__doc__)
dok_array.__doc__ = _matrix_doc_to_array(dok_matrix.__doc__)
lil_array.__doc__ = _matrix_doc_to_array(lil_matrix.__doc__)
