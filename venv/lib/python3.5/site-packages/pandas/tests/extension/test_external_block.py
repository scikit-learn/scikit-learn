# -*- coding: utf-8 -*-
# pylint: disable=W0102

import numpy as np

import pandas as pd
from pandas.core.internals import (
    BlockManager, SingleBlockManager, NonConsolidatableMixIn, Block)

import pytest


class CustomBlock(NonConsolidatableMixIn, Block):

    _holder = np.ndarray

    def formatting_values(self):
        return np.array(["Val: {}".format(i) for i in self.values])

    def concat_same_type(self, to_concat, placement=None):
        """
        Always concatenate disregarding self.ndim as the values are
        always 1D in this custom Block
        """
        values = np.concatenate([blk.values for blk in to_concat])
        return self.make_block_same_class(
            values, placement=placement or slice(0, len(values), 1))


@pytest.fixture
def df():
    df1 = pd.DataFrame({'a': [1, 2, 3]})
    blocks = df1._data.blocks
    values = np.arange(3, dtype='int64')
    custom_block = CustomBlock(values, placement=slice(1, 2))
    blocks = blocks + (custom_block,)
    block_manager = BlockManager(blocks, [pd.Index(['a', 'b']), df1.index])
    return pd.DataFrame(block_manager)


def test_custom_repr():
    values = np.arange(3, dtype='int64')

    # series
    block = CustomBlock(values, placement=slice(0, 3))

    s = pd.Series(SingleBlockManager(block, pd.RangeIndex(3)))
    assert repr(s) == '0    Val: 0\n1    Val: 1\n2    Val: 2\ndtype: int64'

    # dataframe
    block = CustomBlock(values, placement=slice(0, 1))
    blk_mgr = BlockManager([block], [['col'], range(3)])
    df = pd.DataFrame(blk_mgr)
    assert repr(df) == '      col\n0  Val: 0\n1  Val: 1\n2  Val: 2'


def test_concat_series():
    # GH17728
    values = np.arange(3, dtype='int64')
    block = CustomBlock(values, placement=slice(0, 3))
    s = pd.Series(block, pd.RangeIndex(3), fastpath=True)

    res = pd.concat([s, s])
    assert isinstance(res._data.blocks[0], CustomBlock)


def test_concat_dataframe(df):
    # GH17728
    res = pd.concat([df, df])
    assert isinstance(res._data.blocks[1], CustomBlock)


def test_concat_axis1(df):
    # GH17954
    df2 = pd.DataFrame({'c': [.1, .2, .3]})
    res = pd.concat([df, df2], axis=1)
    assert isinstance(res._data.blocks[1], CustomBlock)
