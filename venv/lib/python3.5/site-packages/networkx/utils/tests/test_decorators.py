import tempfile
import os

from nose.tools import *
from nose import SkipTest

import networkx as nx
from networkx.utils.decorators import open_file, not_implemented_for
from networkx.utils.decorators import nodes_or_number, preserve_random_state, \
    random_state


def test_not_implemented_decorator():
    @not_implemented_for('directed')
    def test1(G):
        pass
    test1(nx.Graph())


@raises(KeyError)
def test_not_implemented_decorator_key():
    @not_implemented_for('foo')
    def test1(G):
        pass
    test1(nx.Graph())


@raises(nx.NetworkXNotImplemented)
def test_not_implemented_decorator_raise():
    @not_implemented_for('graph')
    def test1(G):
        pass
    test1(nx.Graph())


class TestOpenFileDecorator(object):
    def setUp(self):
        self.text = ['Blah... ', 'BLAH ', 'BLAH!!!!']
        self.fobj = tempfile.NamedTemporaryFile('wb+', delete=False)
        self.name = self.fobj.name

    def write(self, path):
        for text in self.text:
            path.write(text.encode('ascii'))

    @open_file(1, 'r')
    def read(self, path):
        return path.readlines()[0]

    @staticmethod
    @open_file(0, 'wb')
    def writer_arg0(path):
        path.write('demo'.encode('ascii'))

    @open_file(1, 'wb+')
    def writer_arg1(self, path):
        self.write(path)

    @open_file(2, 'wb')
    def writer_arg2default(self, x, path=None):
        if path is None:
            with tempfile.NamedTemporaryFile('wb+') as fh:
                self.write(fh)
        else:
            self.write(path)

    @open_file(4, 'wb')
    def writer_arg4default(self, x, y, other='hello', path=None, **kwargs):
        if path is None:
            with tempfile.NamedTemporaryFile('wb+') as fh:
                self.write(fh)
        else:
            self.write(path)

    @open_file('path', 'wb')
    def writer_kwarg(self, **kwargs):
        path = kwargs.get('path', None)
        if path is None:
            with tempfile.NamedTemporaryFile('wb+') as fh:
                self.write(fh)
        else:
            self.write(path)

    def test_writer_arg0_str(self):
        self.writer_arg0(self.name)

    def test_writer_arg0_fobj(self):
        self.writer_arg0(self.fobj)

    def test_writer_arg1_str(self):
        self.writer_arg1(self.name)
        assert_equal(self.read(self.name), ''.join(self.text))

    def test_writer_arg1_fobj(self):
        self.writer_arg1(self.fobj)
        assert_false(self.fobj.closed)
        self.fobj.close()
        assert_equal(self.read(self.name), ''.join(self.text))

    def test_writer_arg2default_str(self):
        self.writer_arg2default(0, path=None)
        self.writer_arg2default(0, path=self.name)
        assert_equal(self.read(self.name), ''.join(self.text))

    def test_writer_arg2default_fobj(self):
        self.writer_arg2default(0, path=self.fobj)
        assert_false(self.fobj.closed)
        self.fobj.close()
        assert_equal(self.read(self.name), ''.join(self.text))

    def test_writer_arg2default_fobj(self):
        self.writer_arg2default(0, path=None)

    def test_writer_arg4default_fobj(self):
        self.writer_arg4default(0, 1, dog='dog', other='other')
        self.writer_arg4default(0, 1, dog='dog', other='other', path=self.name)
        assert_equal(self.read(self.name), ''.join(self.text))

    def test_writer_kwarg_str(self):
        self.writer_kwarg(path=self.name)
        assert_equal(self.read(self.name), ''.join(self.text))

    def test_writer_kwarg_fobj(self):
        self.writer_kwarg(path=self.fobj)
        self.fobj.close()
        assert_equal(self.read(self.name), ''.join(self.text))

    def test_writer_kwarg_fobj(self):
        self.writer_kwarg(path=None)

    def tearDown(self):
        self.fobj.close()
        os.unlink(self.name)


@preserve_random_state
def test_preserve_random_state():
    try:
        import numpy.random
        r = numpy.random.random()
    except ImportError:
        return
    assert(abs(r - 0.61879477158568) < 1e-16)


class TestRandomState(object):
    @classmethod
    def setUp(cls):
        global np
        try:
            import numpy as np
        except ImportError:
            raise SkipTest('NumPy not available.')

    @random_state(1)
    def instantiate_random_state(self, random_state):
        assert_true(isinstance(random_state, np.random.RandomState))
        return random_state

    def test_random_state_None(self):
        self.instantiate_random_state(random_state=None)

    def test_random_state_np_random(self):
        self.instantiate_random_state(random_state=np.random)

    def test_random_state_int(self):
        seed = 1
        random_state = self.instantiate_random_state(random_state=seed)
        assert_true(np.all((np.random.RandomState(seed).rand(10),
                            random_state.rand(10))))

    def test_random_state_np_random_RandomState(self):
        seed = 1
        rng = np.random.RandomState(seed)
        random_state = self.instantiate_random_state(random_state=rng)
        assert_true(np.all((np.random.RandomState(seed).rand(10),
                            random_state.rand(10))))


@raises(nx.NetworkXError)
def test_string_arg_index():
    @random_state('a')
    def make_random_state(rs):
        pass
    rstate = make_random_state(1)


@raises(nx.NetworkXError)
def test_invalid_arg_index():
    @random_state(2)
    def make_random_state(rs):
        pass
    rstate = make_random_state(1)
