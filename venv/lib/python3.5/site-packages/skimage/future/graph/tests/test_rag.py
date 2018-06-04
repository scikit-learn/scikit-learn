import numpy as np
from skimage.future import graph
from skimage._shared.version_requirements import is_installed
from skimage import segmentation
from skimage._shared import testing


def max_edge(g, src, dst, n):
    default = {'weight': -np.inf}
    w1 = g[n].get(src, default)['weight']
    w2 = g[n].get(dst, default)['weight']
    return {'weight': max(w1, w2)}


@testing.skipif(not is_installed('networkx'),
                reason="networkx not installed")
def test_rag_merge():
    g = graph.rag.RAG()

    for i in range(5):
        g.add_node(i, {'labels': [i]})

    g.add_edge(0, 1, {'weight': 10})
    g.add_edge(1, 2, {'weight': 20})
    g.add_edge(2, 3, {'weight': 30})
    g.add_edge(3, 0, {'weight': 40})
    g.add_edge(0, 2, {'weight': 50})
    g.add_edge(3, 4, {'weight': 60})

    gc = g.copy()

    # We merge nodes and ensure that the minimum weight is chosen
    # when there is a conflict.
    g.merge_nodes(0, 2)
    assert g.adj[1][2]['weight'] == 10
    assert g.adj[2][3]['weight'] == 30

    # We specify `max_edge` as `weight_func` as ensure that maximum
    # weight is chosen in case on conflict
    gc.merge_nodes(0, 2, weight_func=max_edge)
    assert gc.adj[1][2]['weight'] == 20
    assert gc.adj[2][3]['weight'] == 40

    g.merge_nodes(1, 4)
    g.merge_nodes(2, 3)
    n = g.merge_nodes(3, 4, in_place=False)
    assert sorted(g.node[n]['labels']) == list(range(5))
    assert list(g.edges()) == []


@testing.skipif(not is_installed('networkx'),
                reason="networkx not installed")
def test_threshold_cut():

    img = np.zeros((100, 100, 3), dtype='uint8')
    img[:50, :50] = 255, 255, 255
    img[:50, 50:] = 254, 254, 254
    img[50:, :50] = 2, 2, 2
    img[50:, 50:] = 1, 1, 1

    labels = np.zeros((100, 100), dtype='uint8')
    labels[:50, :50] = 0
    labels[:50, 50:] = 1
    labels[50:, :50] = 2
    labels[50:, 50:] = 3

    rag = graph.rag_mean_color(img, labels)
    new_labels = graph.cut_threshold(labels, rag, 10, in_place=False)
    # Two labels
    assert new_labels.max() == 1

    new_labels = graph.cut_threshold(labels, rag, 10)
    # Two labels
    assert new_labels.max() == 1


@testing.skipif(not is_installed('networkx'),
                reason="networkx not installed")
def test_cut_normalized():

    img = np.zeros((100, 100, 3), dtype='uint8')
    img[:50, :50] = 255, 255, 255
    img[:50, 50:] = 254, 254, 254
    img[50:, :50] = 2, 2, 2
    img[50:, 50:] = 1, 1, 1

    labels = np.zeros((100, 100), dtype='uint8')
    labels[:50, :50] = 0
    labels[:50, 50:] = 1
    labels[50:, :50] = 2
    labels[50:, 50:] = 3

    rag = graph.rag_mean_color(img, labels, mode='similarity')

    new_labels = graph.cut_normalized(labels, rag, in_place=False)
    new_labels, _, _ = segmentation.relabel_sequential(new_labels)
    # Two labels
    assert new_labels.max() == 1

    new_labels = graph.cut_normalized(labels, rag)
    new_labels, _, _ = segmentation.relabel_sequential(new_labels)
    assert new_labels.max() == 1


@testing.skipif(not is_installed('networkx'),
                reason="networkx not installed")
def test_rag_error():
    img = np.zeros((10, 10, 3), dtype='uint8')
    labels = np.zeros((10, 10), dtype='uint8')
    labels[:5, :] = 0
    labels[5:, :] = 1
    with testing.raises(ValueError):
        graph.rag_mean_color(img, labels,
                             2, 'non existant mode')


def _weight_mean_color(graph, src, dst, n):
    diff = graph.node[dst]['mean color'] - graph.node[n]['mean color']
    diff = np.linalg.norm(diff)
    return {'weight': diff}


def _pre_merge_mean_color(graph, src, dst):
    graph.node[dst]['total color'] += graph.node[src]['total color']
    graph.node[dst]['pixel count'] += graph.node[src]['pixel count']
    graph.node[dst]['mean color'] = (graph.node[dst]['total color'] /
                                     graph.node[dst]['pixel count'])


def merge_hierarchical_mean_color(labels, rag, thresh, rag_copy=True,
                                  in_place_merge=False):
    return graph.merge_hierarchical(labels, rag, thresh, rag_copy,
                                    in_place_merge, _pre_merge_mean_color,
                                    _weight_mean_color)


@testing.skipif(not is_installed('networkx'),
                reason="networkx not installed")
def test_rag_hierarchical():
    img = np.zeros((8, 8, 3), dtype='uint8')
    labels = np.zeros((8, 8), dtype='uint8')

    img[:, :, :] = 31
    labels[:, :] = 1

    img[0:4, 0:4, :] = 10, 10, 10
    labels[0:4, 0:4] = 2

    img[4:, 0:4, :] = 20, 20, 20
    labels[4:, 0:4] = 3

    g = graph.rag_mean_color(img, labels)
    g2 = g.copy()
    thresh = 20  # more than 11*sqrt(3) but less than

    result = merge_hierarchical_mean_color(labels, g, thresh)
    assert(np.all(result[:, :4] == result[0, 0]))
    assert(np.all(result[:, 4:] == result[-1, -1]))

    result = merge_hierarchical_mean_color(labels, g2, thresh,
                                           in_place_merge=True)
    assert(np.all(result[:, :4] == result[0, 0]))
    assert(np.all(result[:, 4:] == result[-1, -1]))

    result = graph.cut_threshold(labels, g, thresh)
    assert np.all(result == result[0, 0])


@testing.skipif(not is_installed('networkx'),
                reason="networkx not installed")
def test_ncut_stable_subgraph():
    """ Test to catch an error thrown when subgraph has all equal edges. """

    img = np.zeros((100, 100, 3), dtype='uint8')

    labels = np.zeros((100, 100), dtype='uint8')
    labels[:50, :50] = 1
    labels[:50, 50:] = 2

    rag = graph.rag_mean_color(img, labels, mode='similarity')
    new_labels = graph.cut_normalized(labels, rag, in_place=False)
    new_labels, _, _ = segmentation.relabel_sequential(new_labels)

    assert new_labels.max() == 0


def test_generic_rag_2d():
    labels = np.array([[1, 2], [3, 4]], dtype=np.uint8)
    g = graph.RAG(labels)
    assert g.has_edge(1, 2) and g.has_edge(2, 4) and not g.has_edge(1, 4)
    h = graph.RAG(labels, connectivity=2)
    assert h.has_edge(1, 2) and h.has_edge(1, 4) and h.has_edge(2, 3)


def test_generic_rag_3d():
    labels = np.arange(8, dtype=np.uint8).reshape((2, 2, 2))
    g = graph.RAG(labels)
    assert g.has_edge(0, 1) and g.has_edge(1, 3) and not g.has_edge(0, 3)
    h = graph.RAG(labels, connectivity=2)
    assert h.has_edge(0, 1) and h.has_edge(0, 3) and not h.has_edge(0, 7)
    k = graph.RAG(labels, connectivity=3)
    assert k.has_edge(0, 1) and k.has_edge(1, 2) and k.has_edge(2, 5)


def test_rag_boundary():
    labels = np.zeros((16, 16), dtype='uint8')
    edge_map = np.zeros_like(labels, dtype=float)

    edge_map[8, :] = 0.5
    edge_map[:, 8] = 1.0

    labels[:8, :8] = 1
    labels[:8, 8:] = 2
    labels[8:, :8] = 3
    labels[8:, 8:] = 4

    g = graph.rag_boundary(labels, edge_map, connectivity=1)
    assert set(g.nodes()) == set([1, 2, 3, 4])
    assert set(g.edges()) == set([(1, 2), (1, 3), (2, 4), (3, 4)])
    assert g[1][3]['weight'] == 0.25
    assert g[2][4]['weight'] == 0.34375
    assert g[1][3]['count'] == 16
