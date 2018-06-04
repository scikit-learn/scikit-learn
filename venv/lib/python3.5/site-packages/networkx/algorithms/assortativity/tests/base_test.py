import networkx as nx


class BaseTestAttributeMixing(object):

    def setUp(self):
        G = nx.Graph()
        G.add_nodes_from([0, 1], fish='one')
        G.add_nodes_from([2, 3], fish='two')
        G.add_nodes_from([4], fish='red')
        G.add_nodes_from([5], fish='blue')
        G.add_edges_from([(0, 1), (2, 3), (0, 4), (2, 5)])
        self.G = G

        D = nx.DiGraph()
        D.add_nodes_from([0, 1], fish='one')
        D.add_nodes_from([2, 3], fish='two')
        D.add_nodes_from([4], fish='red')
        D.add_nodes_from([5], fish='blue')
        D.add_edges_from([(0, 1), (2, 3), (0, 4), (2, 5)])
        self.D = D

        M = nx.MultiGraph()
        M.add_nodes_from([0, 1], fish='one')
        M.add_nodes_from([2, 3], fish='two')
        M.add_nodes_from([4], fish='red')
        M.add_nodes_from([5], fish='blue')
        M.add_edges_from([(0, 1), (0, 1), (2, 3)])
        self.M = M

        S = nx.Graph()
        S.add_nodes_from([0, 1], fish='one')
        S.add_nodes_from([2, 3], fish='two')
        S.add_nodes_from([4], fish='red')
        S.add_nodes_from([5], fish='blue')
        S.add_edge(0, 0)
        S.add_edge(2, 2)
        self.S = S


class BaseTestDegreeMixing(object):

    def setUp(self):
        self.P4 = nx.path_graph(4)
        self.D = nx.DiGraph()
        self.D.add_edges_from([(0, 2), (0, 3), (1, 3), (2, 3)])
        self.M = nx.MultiGraph()
        nx.add_path(self.M, range(4))
        self.M.add_edge(0, 1)
        self.S = nx.Graph()
        self.S.add_edges_from([(0, 0), (1, 1)])
