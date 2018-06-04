#!/usr/bin/env python
from nose.tools import *
from nose import SkipTest
import networkx as nx
from base_test import BaseTestAttributeMixing, BaseTestDegreeMixing


class TestDegreeMixingDict(BaseTestDegreeMixing):

    def test_degree_mixing_dict_undirected(self):
        d = nx.degree_mixing_dict(self.P4)
        d_result = {1: {2: 2},
                    2: {1: 2, 2: 2},
                    }
        assert_equal(d, d_result)

    def test_degree_mixing_dict_undirected_normalized(self):
        d = nx.degree_mixing_dict(self.P4, normalized=True)
        d_result = {1: {2: 1.0 / 3},
                    2: {1: 1.0 / 3, 2: 1.0 / 3},
                    }
        assert_equal(d, d_result)

    def test_degree_mixing_dict_directed(self):
        d = nx.degree_mixing_dict(self.D)
        print(d)
        d_result = {1: {3: 2},
                    2: {1: 1, 3: 1},
                    3: {}
                    }
        assert_equal(d, d_result)

    def test_degree_mixing_dict_multigraph(self):
        d = nx.degree_mixing_dict(self.M)
        d_result = {1: {2: 1},
                    2: {1: 1, 3: 3},
                    3: {2: 3}
                    }
        assert_equal(d, d_result)


class TestDegreeMixingMatrix(BaseTestDegreeMixing):

    @classmethod
    def setupClass(cls):
        global np
        global npt
        try:
            import numpy as np
            import numpy.testing as npt

        except ImportError:
            raise SkipTest('NumPy not available.')

    def test_degree_mixing_matrix_undirected(self):
        a_result = np.array([[0, 0, 0],
                             [0, 0, 2],
                             [0, 2, 2]]
                            )
        a = nx.degree_mixing_matrix(self.P4, normalized=False)
        npt.assert_equal(a, a_result)
        a = nx.degree_mixing_matrix(self.P4)
        npt.assert_equal(a, a_result / float(a_result.sum()))

    def test_degree_mixing_matrix_directed(self):
        a_result = np.array([[0, 0, 0, 0],
                             [0, 0, 0, 2],
                             [0, 1, 0, 1],
                             [0, 0, 0, 0]]
                            )
        a = nx.degree_mixing_matrix(self.D, normalized=False)
        npt.assert_equal(a, a_result)
        a = nx.degree_mixing_matrix(self.D)
        npt.assert_equal(a, a_result / float(a_result.sum()))

    def test_degree_mixing_matrix_multigraph(self):
        a_result = np.array([[0, 0, 0, 0],
                             [0, 0, 1, 0],
                             [0, 1, 0, 3],
                             [0, 0, 3, 0]]
                            )
        a = nx.degree_mixing_matrix(self.M, normalized=False)
        npt.assert_equal(a, a_result)
        a = nx.degree_mixing_matrix(self.M)
        npt.assert_equal(a, a_result / float(a_result.sum()))

    def test_degree_mixing_matrix_selfloop(self):
        a_result = np.array([[0, 0, 0],
                             [0, 0, 0],
                             [0, 0, 2]]
                            )
        a = nx.degree_mixing_matrix(self.S, normalized=False)
        npt.assert_equal(a, a_result)
        a = nx.degree_mixing_matrix(self.S)
        npt.assert_equal(a, a_result / float(a_result.sum()))


class TestAttributeMixingDict(BaseTestAttributeMixing):

    def test_attribute_mixing_dict_undirected(self):
        d = nx.attribute_mixing_dict(self.G, 'fish')
        d_result = {'one': {'one': 2, 'red': 1},
                    'two': {'two': 2, 'blue': 1},
                    'red': {'one': 1},
                    'blue': {'two': 1}
                    }
        assert_equal(d, d_result)

    def test_attribute_mixing_dict_directed(self):
        d = nx.attribute_mixing_dict(self.D, 'fish')
        d_result = {'one': {'one': 1, 'red': 1},
                    'two': {'two': 1, 'blue': 1},
                    'red': {},
                    'blue': {}
                    }
        assert_equal(d, d_result)

    def test_attribute_mixing_dict_multigraph(self):
        d = nx.attribute_mixing_dict(self.M, 'fish')
        d_result = {'one': {'one': 4},
                    'two': {'two': 2},
                    }
        assert_equal(d, d_result)


class TestAttributeMixingMatrix(BaseTestAttributeMixing):
    @classmethod
    def setupClass(cls):
        global np
        global npt
        try:
            import numpy as np
            import numpy.testing as npt

        except ImportError:
            raise SkipTest('NumPy not available.')

    def test_attribute_mixing_matrix_undirected(self):
        mapping = {'one': 0, 'two': 1, 'red': 2, 'blue': 3}
        a_result = np.array([[2, 0, 1, 0],
                             [0, 2, 0, 1],
                             [1, 0, 0, 0],
                             [0, 1, 0, 0]]
                            )
        a = nx.attribute_mixing_matrix(self.G, 'fish',
                                       mapping=mapping,
                                       normalized=False)
        npt.assert_equal(a, a_result)
        a = nx.attribute_mixing_matrix(self.G, 'fish',
                                       mapping=mapping)
        npt.assert_equal(a, a_result / float(a_result.sum()))

    def test_attribute_mixing_matrix_directed(self):
        mapping = {'one': 0, 'two': 1, 'red': 2, 'blue': 3}
        a_result = np.array([[1, 0, 1, 0],
                             [0, 1, 0, 1],
                             [0, 0, 0, 0],
                             [0, 0, 0, 0]]
                            )
        a = nx.attribute_mixing_matrix(self.D, 'fish',
                                       mapping=mapping,
                                       normalized=False)
        npt.assert_equal(a, a_result)
        a = nx.attribute_mixing_matrix(self.D, 'fish',
                                       mapping=mapping)
        npt.assert_equal(a, a_result / float(a_result.sum()))

    def test_attribute_mixing_matrix_multigraph(self):
        mapping = {'one': 0, 'two': 1, 'red': 2, 'blue': 3}
        a_result = np.array([[4, 0, 0, 0],
                             [0, 2, 0, 0],
                             [0, 0, 0, 0],
                             [0, 0, 0, 0]]
                            )
        a = nx.attribute_mixing_matrix(self.M, 'fish',
                                       mapping=mapping,
                                       normalized=False)
        npt.assert_equal(a, a_result)
        a = nx.attribute_mixing_matrix(self.M, 'fish',
                                       mapping=mapping)
        npt.assert_equal(a, a_result / float(a_result.sum()))
