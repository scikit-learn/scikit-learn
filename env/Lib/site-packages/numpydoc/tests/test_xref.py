import pytest

from numpydoc.xref import DEFAULT_LINKS, make_xref

# Use the default numpydoc link mapping
xref_aliases = DEFAULT_LINKS


# Comes mainly from numpy
data = r"""
(...) array_like, float, optional
(...) :term:`numpy:array_like`, :class:`python:float`, optional

(2,) ndarray
(2,) :obj:`ndarray <numpy.ndarray>`

(...,M,N) array_like
(...,M,N) :term:`numpy:array_like`

(..., M, N) array_like
(..., :obj:`M`, :obj:`N`) :term:`numpy:array_like`

(float, float), optional
(:class:`python:float`, :class:`python:float`), optional

1-D array or sequence
1-D :obj:`array <numpy.ndarray>` or :term:`python:sequence`

array of str or unicode-like
:obj:`array <numpy.ndarray>` of :class:`python:str` or unicode-like

array_like of float
:term:`numpy:array_like` of :class:`python:float`

bool or callable
:ref:`bool <python:bltin-boolean-values>` or :func:`python:callable`

int in [0, 255]
:class:`python:int` in [0, 255]

int or None, optional
:class:`python:int` or :data:`python:None`, optional

list of str or array_like
:class:`python:list` of :class:`python:str` or :term:`numpy:array_like`

sequence of array_like
:term:`python:sequence` of :term:`numpy:array_like`

str or pathlib.Path
:class:`python:str` or :obj:`pathlib.Path`

{'', string}, optional
{'', :class:`python:str`}, optional

{'C', 'F', 'A', or 'K'}, optional
{'C', 'F', 'A', or 'K'}, optional

{'linear', 'lower', 'higher', 'midpoint', 'nearest'}
{'linear', 'lower', 'higher', 'midpoint', 'nearest'}

{False, True, 'greedy', 'optimal'}
{:data:`python:False`, :data:`python:True`, 'greedy', 'optimal'}

{{'begin', 1}, {'end', 0}}, {string, int}
{{'begin', 1}, {'end', 0}}, {:class:`python:str`, :class:`python:int`}

callable f'(x,*args)
:func:`python:callable` f'(x,*args)

callable ``fhess(x, *args)``, optional
:func:`python:callable` ``fhess(x, *args)``, optional

spmatrix (format: ``csr``, ``bsr``, ``dia`` or coo``)
:obj:`spmatrix` (format: ``csr``, ``bsr``, ``dia`` or coo``)

:ref:`strftime <strftime-strptime-behavior>`
:ref:`strftime <strftime-strptime-behavior>`

callable or :ref:`strftime <strftime-strptime-behavior>`
:func:`python:callable` or :ref:`strftime <strftime-strptime-behavior>`

callable or :ref:`strftime behavior <strftime-strptime-behavior>`
:func:`python:callable` or :ref:`strftime behavior <strftime-strptime-behavior>`

list(int)
:class:`python:list`\(:class:`python:int`)

list[int]
:class:`python:list`\[:class:`python:int`]

dict(str, int)
:class:`python:dict`\(:class:`python:str`, :class:`python:int`)

dict[str,  int]
:class:`python:dict`\[:class:`python:str`,  :class:`python:int`]

tuple(float, float)
:class:`python:tuple`\(:class:`python:float`, :class:`python:float`)

dict[tuple(str, str), int]
:class:`python:dict`\[:class:`python:tuple`\(:class:`python:str`, :class:`python:str`), :class:`python:int`]
"""

data_ignore_obj = r"""
(...) array_like, float, optional
(...) :term:`numpy:array_like`, :class:`python:float`, optional

(2,) ndarray
(2,) :obj:`ndarray <numpy.ndarray>`

(...,M,N) array_like
(...,M,N) :term:`numpy:array_like`

(..., M, N) array_like
(..., M, N) :term:`numpy:array_like`

(float, float), optional
(:class:`python:float`, :class:`python:float`), optional

1-D array or sequence
1-D :obj:`array <numpy.ndarray>` or :term:`python:sequence`

array of str or unicode-like
:obj:`array <numpy.ndarray>` of :class:`python:str` or unicode-like

array_like of float
:term:`numpy:array_like` of :class:`python:float`

bool or callable
:ref:`bool <python:bltin-boolean-values>` or :func:`python:callable`

int in [0, 255]
:class:`python:int` in [0, 255]

int or None, optional
:class:`python:int` or :data:`python:None`, optional

list of str or array_like
:class:`python:list` of :class:`python:str` or :term:`numpy:array_like`

sequence of array_like
:term:`python:sequence` of :term:`numpy:array_like`

str or pathlib.Path
:class:`python:str` or pathlib.Path

{'', string}, optional
{'', :class:`python:str`}, optional

{'C', 'F', 'A', or 'K'}, optional
{'C', 'F', 'A', or 'K'}, optional

{'linear', 'lower', 'higher', 'midpoint', 'nearest'}
{'linear', 'lower', 'higher', 'midpoint', 'nearest'}

{False, True, 'greedy', 'optimal'}
{:data:`python:False`, :data:`python:True`, 'greedy', 'optimal'}

{{'begin', 1}, {'end', 0}}, {string, int}
{{'begin', 1}, {'end', 0}}, {:class:`python:str`, :class:`python:int`}

callable f'(x,*args)
:func:`python:callable` f'(x,*args)

callable ``fhess(x, *args)``, optional
:func:`python:callable` ``fhess(x, *args)``, optional

spmatrix (format: ``csr``, ``bsr``, ``dia`` or coo``)
spmatrix (format: ``csr``, ``bsr``, ``dia`` or coo``)

:ref:`strftime <strftime-strptime-behavior>`
:ref:`strftime <strftime-strptime-behavior>`

callable or :ref:`strftime <strftime-strptime-behavior>`
:func:`python:callable` or :ref:`strftime <strftime-strptime-behavior>`

callable or :ref:`strftime behavior <strftime-strptime-behavior>`
:func:`python:callable` or :ref:`strftime behavior <strftime-strptime-behavior>`

list(int)
:class:`python:list`\(:class:`python:int`)

list[int]
:class:`python:list`\[:class:`python:int`]

dict(str, int)
:class:`python:dict`\(:class:`python:str`, :class:`python:int`)

dict[str,  int]
:class:`python:dict`\[:class:`python:str`,  :class:`python:int`]

tuple(float, float)
:class:`python:tuple`\(:class:`python:float`, :class:`python:float`)

dict[tuple(str, str), int]
:class:`python:dict`\[:class:`python:tuple`\(:class:`python:str`, :class:`python:str`), :class:`python:int`]
"""

xref_ignore = {"or", "in", "of", "default", "optional"}


@pytest.mark.parametrize(
    ("param_type", "expected_result"),
    [tuple(s.split("\n")) for s in data.strip().split("\n\n")],
)
def test_make_xref(param_type, expected_result):
    assert make_xref(param_type, xref_aliases, xref_ignore) == expected_result


@pytest.mark.parametrize(
    ("param_type", "expected_result"),
    [tuple(s.split("\n")) for s in data_ignore_obj.strip().split("\n\n")],
)
def test_make_xref_ignore_unknown(param_type, expected_result):
    assert make_xref(param_type, xref_aliases, xref_ignore="all") == expected_result


def test_xref_ignore_is_all():
    with pytest.raises(TypeError, match="must be a set or 'all'"):
        make_xref("array_like", xref_aliases, xref_ignore="foo")
