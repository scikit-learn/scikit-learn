#!/usr/bin/env python
from nose.tools import *
from networkx.utils import powerlaw_sequence,\
    zipf_rv, random_weighted_sample,\
    weighted_choice
import networkx.utils


def test_degree_sequences():
    seq = powerlaw_sequence(10)
    assert_equal(len(seq), 10)


def test_zipf_rv():
    r = zipf_rv(2.3)
    assert_true(type(r), int)
    assert_raises(ValueError, zipf_rv, 0.5)
    assert_raises(ValueError, zipf_rv, 2, xmin=0)


def test_random_weighted_sample():
    mapping = {'a': 10, 'b': 20}
    s = random_weighted_sample(mapping, 2)
    assert_equal(sorted(s), sorted(mapping.keys()))
    assert_raises(ValueError, random_weighted_sample, mapping, 3)


def test_random_weighted_choice():
    mapping = {'a': 10, 'b': 0}
    c = weighted_choice(mapping)
    assert_equal(c, 'a')
