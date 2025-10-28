"""Numpydoc test module.

.. currentmodule:: numpydoc_test_module

.. autosummary::
   :toctree: generated/

   MyClass
   my_function

Reference [1]_

References
----------
.. [1] https://numpydoc.readthedocs.io
"""

__all__ = ["MyClass", "my_function"]


class MyClass:
    """A class.

    Reference [2]_

    Parameters
    ----------
    *args : iterable
        Arguments.
    **kwargs : dict
        Keyword arguments.

    References
    ----------
    .. [2] https://numpydoc.readthedocs.io
    """

    def example(self, x):
        """Example method."""


def my_function(*args, **kwargs):
    """Return None.

    See [3]_.

    Parameters
    ----------
    *args : iterable
        Arguments.
    **kwargs : dict
        Keyword arguments.

    Returns
    -------
    out : None
        The output.

    References
    ----------
    .. [3] https://numpydoc.readthedocs.io
    """
    return
