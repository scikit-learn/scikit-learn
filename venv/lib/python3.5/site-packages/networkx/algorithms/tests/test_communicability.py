from collections import defaultdict
from nose.tools import *
from nose import SkipTest
import networkx as nx
from networkx.algorithms.communicability_alg import *


class TestCommunicability:
    @classmethod
    def setupClass(cls):
        global numpy
        global scipy
        try:
            import numpy
        except ImportError:
            raise SkipTest('NumPy not available.')
        try:
            import scipy
        except ImportError:
            raise SkipTest('SciPy not available.')

    def test_communicability(self):
        answer = {0: {0: 1.5430806348152435,
                      1: 1.1752011936438012
                      },
                  1: {0: 1.1752011936438012,
                      1: 1.5430806348152435
                      }
                  }
#        answer={(0, 0): 1.5430806348152435,
#                (0, 1): 1.1752011936438012,
#                (1, 0): 1.1752011936438012,
#                (1, 1): 1.5430806348152435}

        result = communicability(nx.path_graph(2))
        for k1, val in result.items():
            for k2 in val:
                assert_almost_equal(answer[k1][k2], result[k1][k2], places=7)

    def test_communicability2(self):

        answer_orig = {('1', '1'): 1.6445956054135658,
                       ('1', 'Albert'): 0.7430186221096251,
                       ('1', 'Aric'): 0.7430186221096251,
                       ('1', 'Dan'): 1.6208126320442937,
                       ('1', 'Franck'): 0.42639707170035257,
                       ('Albert', '1'): 0.7430186221096251,
                       ('Albert', 'Albert'): 2.4368257358712189,
                       ('Albert', 'Aric'): 1.4368257358712191,
                       ('Albert', 'Dan'): 2.0472097037446453,
                       ('Albert', 'Franck'): 1.8340111678944691,
                       ('Aric', '1'): 0.7430186221096251,
                       ('Aric', 'Albert'): 1.4368257358712191,
                       ('Aric', 'Aric'): 2.4368257358712193,
                       ('Aric', 'Dan'): 2.0472097037446457,
                       ('Aric', 'Franck'): 1.8340111678944691,
                       ('Dan', '1'): 1.6208126320442937,
                       ('Dan', 'Albert'): 2.0472097037446453,
                       ('Dan', 'Aric'): 2.0472097037446457,
                       ('Dan', 'Dan'): 3.1306328496328168,
                       ('Dan', 'Franck'): 1.4860372442192515,
                       ('Franck', '1'): 0.42639707170035257,
                       ('Franck', 'Albert'): 1.8340111678944691,
                       ('Franck', 'Aric'): 1.8340111678944691,
                       ('Franck', 'Dan'): 1.4860372442192515,
                       ('Franck', 'Franck'): 2.3876142275231915}

        answer = defaultdict(dict)
        for (k1, k2), v in answer_orig.items():
            answer[k1][k2] = v

        G1 = nx.Graph([('Franck', 'Aric'), ('Aric', 'Dan'), ('Dan', 'Albert'),
                       ('Albert', 'Franck'), ('Dan', '1'), ('Franck', 'Albert')])

        result = communicability(G1)
        for k1, val in result.items():
            for k2 in val:
                assert_almost_equal(answer[k1][k2], result[k1][k2], places=7)

        result = communicability_exp(G1)
        for k1, val in result.items():
            for k2 in val:
                assert_almost_equal(answer[k1][k2], result[k1][k2], places=7)
