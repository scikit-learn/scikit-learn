import os
from functools import partial
import re
from operator import add, neg
import sys
import pytest

optimize2 = (sys.flags.optimize == 2)
if not optimize2:
    pytest.importorskip("graphviz")
    from dask.dot import dot_graph, task_label, label, to_graphviz
else:
    pytestmark = pytest.mark.skipif(True,
                                    reason="graphviz exception with Python -OO flag")

from dask import delayed
from dask.utils import ensure_not_exists
from IPython.display import Image, SVG


# Since graphviz doesn't store a graph, we need to parse the output
label_re = re.compile('.*\[label=(.*?) shape=(.*?)\]')


def get_label(line):
    m = label_re.match(line)
    if m:
        return m.group(1)


def get_shape(line):
    m = label_re.match(line)
    if m:
        return m.group(2)


dsk = {'a': 1,
       'b': 2,
       'c': (neg, 'a'),
       'd': (neg, 'b'),
       'e': (add, 'c', 'd'),
       'f': (sum, ['a', 'e'])}


def test_task_label():
    assert task_label((partial(add, 1), 1)) == 'add'
    assert task_label((add, 1)) == 'add'
    assert task_label((add, (add, 1, 2))) == 'add(...)'


def test_label():
    assert label('x') == 'x'
    assert label('elemwise-ffcd9aa2231d466b5aa91e8bfa9e9487') == 'elemwise-#'

    cache = {}
    result = label('elemwise-ffcd9aa2231d466b5aa91e8bfa9e9487', cache=cache)
    assert result == 'elemwise-#0'
    # cached
    result = label('elemwise-ffcd9aa2231d466b5aa91e8bfa9e9487', cache=cache)
    assert result == 'elemwise-#0'
    assert len(cache) == 1

    result = label('elemwise-e890b510984f344edea9a5e5fe05c0db', cache=cache)
    assert result == 'elemwise-#1'
    assert len(cache) == 2

    result = label('elemwise-ffcd9aa2231d466b5aa91e8bfa9e9487', cache=cache)
    assert result == 'elemwise-#0'
    assert len(cache) == 2

    assert label('x', cache=cache) == 'x'
    assert len(cache) == 2


def test_to_graphviz():
    g = to_graphviz(dsk)
    labels = list(filter(None, map(get_label, g.body)))
    assert len(labels) == 10        # 10 nodes total
    funcs = set(('add', 'sum', 'neg'))
    assert set(labels).difference(dsk) == funcs
    assert set(labels).difference(funcs) == set(dsk)
    shapes = list(filter(None, map(get_shape, g.body)))
    assert set(shapes) == set(('box', 'circle'))


def test_to_graphviz_custom():
    g = to_graphviz(
        dsk,
        data_attributes={'a': {'shape': 'square'}},
        function_attributes={'c': {'label': 'neg_c', 'shape': 'ellipse'}},
    )
    labels = list(filter(None, map(get_label, g.body)))
    funcs = set(('add', 'sum', 'neg', 'neg_c'))
    assert set(labels).difference(dsk) == funcs
    assert set(labels).difference(funcs) == set(dsk)
    shapes = list(filter(None, map(get_shape, g.body)))
    assert set(shapes) == set(('box', 'circle', 'square', 'ellipse'))


def test_to_graphviz_attributes():
    assert to_graphviz(dsk).graph_attr['rankdir'] == 'BT'
    assert to_graphviz(dsk, rankdir='LR').graph_attr['rankdir'] == 'LR'
    assert to_graphviz(dsk, node_attr={'color': 'white'}).node_attr['color'] == 'white'
    assert to_graphviz(dsk, edge_attr={'color': 'white'}).edge_attr['color'] == 'white'


def test_aliases():
    g = to_graphviz({'x': 1, 'y': 'x'})
    labels = list(filter(None, map(get_label, g.body)))
    assert len(labels) == 2
    assert len(g.body) - len(labels) == 1   # Single edge


def test_dot_graph(tmpdir):
    # Use a name that the shell would interpret specially to ensure that we're
    # not vulnerable to shell injection when interacting with `dot`.
    filename = str(tmpdir.join('$(touch should_not_get_created.txt)'))

    # Map from format extension to expected return type.
    result_types = {
        'png': Image,
        'jpeg': Image,
        'dot': type(None),
        'pdf': type(None),
        'svg': SVG,
    }
    for format in result_types:
        target = '.'.join([filename, format])
        ensure_not_exists(target)
        try:
            result = dot_graph(dsk, filename=filename, format=format)

            assert not os.path.exists('should_not_get_created.txt')
            assert os.path.isfile(target)
            assert isinstance(result, result_types[format])
        finally:
            ensure_not_exists(target)


def test_dot_graph_no_filename(tmpdir):
    # Map from format extension to expected return type.
    result_types = {
        'png': Image,
        'jpeg': Image,
        'dot': type(None),
        'pdf': type(None),
        'svg': SVG,
    }
    for format in result_types:
        before = tmpdir.listdir()
        result = dot_graph(dsk, filename=None, format=format)
        # We shouldn't write any files if filename is None.
        after = tmpdir.listdir()
        assert before == after
        assert isinstance(result, result_types[format])


def test_dot_graph_defaults():
    # Test with default args.
    default_name = 'mydask'
    default_format = 'png'
    target = '.'.join([default_name, default_format])

    ensure_not_exists(target)
    try:
        result = dot_graph(dsk)
        assert os.path.isfile(target)
        assert isinstance(result, Image)
    finally:
        ensure_not_exists(target)


def test_filenames_and_formats():
    # Test with a variety of user provided args
    filenames = ['mydaskpdf', 'mydask.pdf', 'mydask.pdf', 'mydaskpdf', 'mydask.pdf.svg']
    formats = ['svg', None, 'svg', None, None]
    targets = ['mydaskpdf.svg', 'mydask.pdf', 'mydask.pdf.svg', 'mydaskpdf.png', 'mydask.pdf.svg']

    result_types = {
        'png': Image,
        'jpeg': Image,
        'dot': type(None),
        'pdf': type(None),
        'svg': SVG,
    }

    for filename, format, target in zip(filenames, formats, targets):
        expected_result_type = result_types[target.split('.')[-1]]
        result = dot_graph(dsk, filename=filename, format=format)
        assert os.path.isfile(target)
        assert isinstance(result, expected_result_type)
        ensure_not_exists(target)


def test_delayed_kwargs_apply():
    def f(x, y=True):
        return x + y

    x = delayed(f)(1, y=2)
    label = task_label(x.dask[x.key])
    assert 'f' in label
    assert 'apply' not in label
