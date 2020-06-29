"""Central place for array creation, to support non-numpy arrays.

This currently leverages NEP18 via np.{empty|zeros|ones}_like to create
non-numpy arrays.
"""

from .fixes import np_version

import numpy as np


def create_like(create, create_like):
    """Generalization of (empty|zeros|ones)_like"""
    name = create.__name__

    def metafunction(prototype, dtype=None, order='C', subok=True, shape=None):
        """Forwards call to numpy.{name}_like or {name}, to be compatible with NEP18.

        Before numpy 1.17, numpy.{name}_like did not take a shape argument.

        When version of numpy < (1, 17), and shape is provided, the call will
        be forwarded to numpy.{name}. If shape is not provided, the call is
        forwarded to numpy.{name}_like.
        """.format(name=name)
        if np_version < (1, 17):
            if shape is not None:
                if dtype is None:
                    if not hasattr(prototype, 'dtype'):
                        raise NotImplementedError('Passed prototype to {name}_'
                                                  'like without a dtype'.
                                                  format(name=name))
                    dtype = prototype.dtype
                if order == 'A':
                    order = 'F' if prototype.flags['F_CONTIGUOUS'] else 'C'
                elif order == 'K':
                    raise NotImplementedError('order=K not implemented')
                return create(shape, dtype=dtype, order=order)
            else:
                return create_like(prototype, dtype=dtype, order=order,
                                   subok=subok)
        else:
            return create_like(prototype, dtype=dtype, order=order,
                               shape=shape)
    return metafunction


empty_like = create_like(np.empty, np.empty_like)
zeros_like = create_like(np.zeros, np.zeros_like)
ones_like = create_like(np.ones, np.ones_like)
