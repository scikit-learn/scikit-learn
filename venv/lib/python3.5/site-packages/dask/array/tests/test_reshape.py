import pytest
import numpy as np
from dask.array.reshape import reshape_rechunk, expand_tuple, contract_tuple


@pytest.mark.parametrize('inshape,outshape,prechunks,inchunks,outchunks', [
    ((4,), (4,), ((2, 2),), ((2, 2),), ((2, 2),)),
    ((4,), (2, 2), ((2, 2),), ((2, 2),), ((1, 1), (2,))),

    ((4,), (4, 1), ((2, 2),), ((2, 2),), ((2, 2), (1,))),
    ((4,), (1, 4), ((2, 2),), ((2, 2),), ((1,), (2, 2))),
    ((1, 4), (4,), ((1,), (2, 2),), ((1,), (2, 2),), ((2, 2),)),
    ((4, 1), (4,), ((2, 2), (1,)), ((2, 2), (1,)), ((2, 2),)),
    ((4, 1, 4), (4, 4), ((2, 2), (1,), (2, 2),), ((2, 2), (1,), (2, 2),), ((2, 2), (2, 2))),
    ((4, 4), (4, 1, 4), ((2, 2), (2, 2),), ((2, 2), (2, 2),), ((2, 2), (1,), (2, 2))),

    ((2, 2), (4,), ((2,), (2,)), ((2,), (2,)), ((4,),)),
    ((2, 2), (4,), ((1, 1), (2,)), ((1, 1,), (2,)), ((2, 2),)),
    ((2, 2), (4,), ((2,), (1, 1)), ((1, 1,), (2,)), ((2, 2),)),

    ((64,), (4, 4, 4), ((8, 8, 8, 8, 8, 8, 8, 8),), ((16, 16, 16, 16),), ((1, 1, 1, 1), (4,), (4,))),
    ((64,), (4, 4, 4), ((32, 32),), ((32, 32),), ((2, 2), (4,), (4,))),
    ((64,), (4, 4, 4), ((16, 48),), ((16, 48),), ((1, 3), (4,), (4,))),
    ((64,), (4, 4, 4), ((20, 44),), ((16, 48),), ((1, 3), (4,), (4,))),
    ((64, 4), (8, 8, 4), ((16, 16, 16, 16), (2, 2)), ((16, 16, 16, 16), (2, 2)), ((2, 2, 2, 2), (8,), (2, 2))),
])
def test_reshape_rechunk(inshape, outshape, prechunks, inchunks, outchunks):
    result_in, result_out = reshape_rechunk(inshape, outshape, prechunks)
    assert result_in == inchunks
    assert result_out == outchunks
    assert np.prod(list(map(len, result_in))) == np.prod(list(map(len, result_out)))


def test_expand_tuple():
    assert expand_tuple((2, 4), 2) == (1, 1, 2, 2)
    assert expand_tuple((2, 4), 3) == (1, 1, 1, 1, 2)
    assert expand_tuple((3, 4), 2) == (1, 2, 2, 2)
    assert expand_tuple((7, 4), 3) == (2, 2, 3, 1, 1, 2)


def test_contract_tuple():
    assert contract_tuple((1, 1, 2, 3, 1), 2) == (2, 2, 2, 2)
    assert contract_tuple((1, 1, 2, 5, 1), 2) == (2, 2, 4, 2)
    assert contract_tuple((2, 4), 2) == (2, 4)
    assert contract_tuple((2, 4), 3) == (6,)
