"""Unit tests for matplotlib drawing functions."""
import os
import itertools
from nose import SkipTest
import networkx as nx


class TestPylab(object):
    @classmethod
    def setupClass(cls):
        global plt
        try:
            import matplotlib as mpl
            mpl.use('PS', warn=False)
            import matplotlib.pyplot as plt
            plt.rcParams['text.usetex'] = False
        except ImportError:
            raise SkipTest('matplotlib not available.')
        except RuntimeError:
            raise SkipTest('matplotlib not available.')

    def setUp(self):
        self.G = nx.barbell_graph(4, 6)

    def test_draw(self):
        try:
            functions = [nx.draw_circular,
                         nx.draw_kamada_kawai,
                         nx.draw_random,
                         nx.draw_spectral,
                         nx.draw_spring,
                         nx.draw_shell]
            options = [{
                'node_color': 'black',
                'node_size': 100,
                'width': 3,
            }]
            for function, option in itertools.product(functions, options):
                function(self.G, **option)
                plt.savefig('test.ps')

        finally:
            try:
                os.unlink('test.ps')
            except OSError:
                pass

    def test_edge_colormap(self):
        colors = range(self.G.number_of_edges())
        nx.draw_spring(self.G, edge_color=colors, width=4,
                       edge_cmap=plt.cm.Blues, with_labels=True)
        plt.show()

    def test_arrows(self):
        nx.draw_spring(self.G.to_directed())
        plt.show()

    def test_edge_colors_and_widths(self):
        nx.draw_random(self.G, edgelist=[(0, 1), (0, 2)], width=[1, 2], edge_colors=['r', 'b'])

    def test_labels_and_colors(self):
        G = nx.cubical_graph()
        pos = nx.spring_layout(G)  # positions for all nodes
        # nodes
        nx.draw_networkx_nodes(G, pos,
                               nodelist=[0, 1, 2, 3],
                               node_color='r',
                               node_size=500,
                               alpha=0.8)
        nx.draw_networkx_nodes(G, pos,
                               nodelist=[4, 5, 6, 7],
                               node_color='b',
                               node_size=500,
                               alpha=0.8)
        # edges
        nx.draw_networkx_edges(G, pos, width=1.0, alpha=0.5)
        nx.draw_networkx_edges(G, pos,
                               edgelist=[(0, 1), (1, 2), (2, 3), (3, 0)],
                               width=8, alpha=0.5, edge_color='r')
        nx.draw_networkx_edges(G, pos,
                               edgelist=[(4, 5), (5, 6), (6, 7), (7, 4)],
                               width=8, alpha=0.5, edge_color='b')
        # some math labels
        labels = {}
        labels[0] = r'$a$'
        labels[1] = r'$b$'
        labels[2] = r'$c$'
        labels[3] = r'$d$'
        labels[4] = r'$\alpha$'
        labels[5] = r'$\beta$'
        labels[6] = r'$\gamma$'
        labels[7] = r'$\delta$'
        nx.draw_networkx_labels(G, pos, labels, font_size=16)
        plt.show()

    def test_axes(self):
        fig, ax = plt.subplots()
        nx.draw(self.G, ax=ax)

    def test_empty_graph(self):
        G = nx.Graph()
        nx.draw(G)

    def test_alpha_iter(self):
        pos = nx.random_layout(self.G)
        # with fewer alpha elements than nodes
        plt.subplot(131)
        nx.draw_networkx_nodes(self.G, pos, alpha=[0.1, 0.2])
        # with equal alpha elements and nodes
        num_nodes = len(self.G.nodes)
        alpha = [x / num_nodes for x in range(num_nodes)]
        colors = range(num_nodes)
        plt.subplot(132)
        nx.draw_networkx_nodes(self.G, pos, node_color=colors, alpha=alpha)
        # with more alpha elements than nodes
        alpha.append(1)
        plt.subplot(133)
        nx.draw_networkx_nodes(self.G, pos, alpha=alpha)
