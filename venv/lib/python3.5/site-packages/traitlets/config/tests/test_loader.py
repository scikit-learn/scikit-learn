# encoding: utf-8
"""Tests for traitlets.config.loader"""

# Copyright (c) IPython Development Team.
# Distributed under the terms of the Modified BSD License.

import copy
import logging
import os
import pickle
import sys
from tempfile import mkstemp
from unittest import TestCase

from pytest import skip

from traitlets.config.loader import (
    Config,
    LazyConfigValue,
    PyFileConfigLoader,
    JSONFileConfigLoader,
    KeyValueConfigLoader,
    ArgParseConfigLoader,
    KVArgParseConfigLoader,
    ConfigError,
)


pyfile = """
c = get_config()
c.a=10
c.b=20
c.Foo.Bar.value=10
c.Foo.Bam.value=list(range(10))
c.D.C.value='hi there'
"""

json1file = """
{
  "version": 1,
  "a": 10,
  "b": 20,
  "Foo": {
    "Bam": {
      "value": [ 0, 1, 2, 3, 4, 5, 6, 7, 8, 9 ]
    },
    "Bar": {
      "value": 10
    }
  },
  "D": {
    "C": {
      "value": "hi there"
    }
  }
}
"""

# should not load
json2file = """
{
  "version": 2
}
"""

import logging
log = logging.getLogger('devnull')
log.setLevel(0)

class TestFileCL(TestCase):

    def _check_conf(self, config):
        self.assertEqual(config.a, 10)
        self.assertEqual(config.b, 20)
        self.assertEqual(config.Foo.Bar.value, 10)
        self.assertEqual(config.Foo.Bam.value, list(range(10)))
        self.assertEqual(config.D.C.value, 'hi there')

    def test_python(self):
        fd, fname = mkstemp('.py')
        f = os.fdopen(fd, 'w')
        f.write(pyfile)
        f.close()
        # Unlink the file
        cl = PyFileConfigLoader(fname, log=log)
        config = cl.load_config()
        self._check_conf(config)

    def test_json(self):
        fd, fname = mkstemp('.json')
        f = os.fdopen(fd, 'w')
        f.write(json1file)
        f.close()
        # Unlink the file
        cl = JSONFileConfigLoader(fname, log=log)
        config = cl.load_config()
        self._check_conf(config)

    def test_context_manager(self):

        fd, fname = mkstemp('.json')
        f = os.fdopen(fd, 'w')
        f.write('{}')
        f.close()

        cl = JSONFileConfigLoader(fname, log=log)

        value = 'context_manager'

        with cl as c:
            c.MyAttr.value = value

        self.assertEqual(cl.config.MyAttr.value, value)

        # check that another loader does see the change
        cl2 = JSONFileConfigLoader(fname, log=log)
        self.assertEqual(cl.config.MyAttr.value, value)

    def test_json_context_bad_write(self):
        fd, fname = mkstemp('.json')
        f = os.fdopen(fd, 'w')
        f.write('{}')
        f.close()
        
        with JSONFileConfigLoader(fname, log=log) as config:
            config.A.b = 1
        
        with self.assertRaises(TypeError):
            with JSONFileConfigLoader(fname, log=log) as config:
                config.A.cant_json = lambda x: x
        
        loader = JSONFileConfigLoader(fname, log=log)
        cfg = loader.load_config()
        assert cfg.A.b == 1
        assert 'cant_json' not in cfg.A

    def test_collision(self):
        a = Config()
        b = Config()
        self.assertEqual(a.collisions(b), {})
        a.A.trait1 = 1
        b.A.trait2 = 2
        self.assertEqual(a.collisions(b), {})
        b.A.trait1 = 1
        self.assertEqual(a.collisions(b), {})
        b.A.trait1 = 0
        self.assertEqual(a.collisions(b), {
            'A': {
                'trait1': "1 ignored, using 0",
            }
        })
        self.assertEqual(b.collisions(a), {
            'A': {
                'trait1': "0 ignored, using 1",
            }
        })
        a.A.trait2 = 3
        self.assertEqual(b.collisions(a), {
            'A': {
                'trait1': "0 ignored, using 1",
                'trait2': "2 ignored, using 3",
            }
        })

    def test_v2raise(self):
        fd, fname = mkstemp('.json')
        f = os.fdopen(fd, 'w')
        f.write(json2file)
        f.close()
        # Unlink the file
        cl = JSONFileConfigLoader(fname, log=log)
        with self.assertRaises(ValueError):
            cl.load_config()


class MyLoader1(ArgParseConfigLoader):
    def _add_arguments(self, aliases=None, flags=None):
        p = self.parser
        p.add_argument('-f', '--foo', dest='Global.foo', type=str)
        p.add_argument('-b', dest='MyClass.bar', type=int)
        p.add_argument('-n', dest='n', action='store_true')
        p.add_argument('Global.bam', type=str)

class MyLoader2(ArgParseConfigLoader):
    def _add_arguments(self, aliases=None, flags=None):
        subparsers = self.parser.add_subparsers(dest='subparser_name')
        subparser1 = subparsers.add_parser('1')
        subparser1.add_argument('-x',dest='Global.x')
        subparser2 = subparsers.add_parser('2')
        subparser2.add_argument('y')

class TestArgParseCL(TestCase):

    def test_basic(self):
        cl = MyLoader1()
        config = cl.load_config('-f hi -b 10 -n wow'.split())
        self.assertEqual(config.Global.foo, 'hi')
        self.assertEqual(config.MyClass.bar, 10)
        self.assertEqual(config.n, True)
        self.assertEqual(config.Global.bam, 'wow')
        config = cl.load_config(['wow'])
        self.assertEqual(list(config.keys()), ['Global'])
        self.assertEqual(list(config.Global.keys()), ['bam'])
        self.assertEqual(config.Global.bam, 'wow')

    def test_add_arguments(self):
        cl = MyLoader2()
        config = cl.load_config('2 frobble'.split())
        self.assertEqual(config.subparser_name, '2')
        self.assertEqual(config.y, 'frobble')
        config = cl.load_config('1 -x frobble'.split())
        self.assertEqual(config.subparser_name, '1')
        self.assertEqual(config.Global.x, 'frobble')

    def test_argv(self):
        cl = MyLoader1(argv='-f hi -b 10 -n wow'.split())
        config = cl.load_config()
        self.assertEqual(config.Global.foo, 'hi')
        self.assertEqual(config.MyClass.bar, 10)
        self.assertEqual(config.n, True)
        self.assertEqual(config.Global.bam, 'wow')


class TestKeyValueCL(TestCase):
    klass = KeyValueConfigLoader

    def test_eval(self):
        cl = self.klass(log=log)
        config = cl.load_config('--Class.str_trait=all --Class.int_trait=5 --Class.list_trait=["hello",5]'.split())
        self.assertEqual(config.Class.str_trait, 'all')
        self.assertEqual(config.Class.int_trait, 5)
        self.assertEqual(config.Class.list_trait, ["hello", 5])

    def test_basic(self):
        cl = self.klass(log=log)
        argv = [ '--' + s[2:] for s in pyfile.split('\n') if s.startswith('c.') ]
        print(argv)
        config = cl.load_config(argv)
        self.assertEqual(config.a, 10)
        self.assertEqual(config.b, 20)
        self.assertEqual(config.Foo.Bar.value, 10)
        # non-literal expressions are not evaluated
        self.assertEqual(config.Foo.Bam.value, 'list(range(10))')
        self.assertEqual(config.D.C.value, 'hi there')
    
    def test_expanduser(self):
        cl = self.klass(log=log)
        argv = ['--a=~/1/2/3', '--b=~', '--c=~/', '--d="~/"']
        config = cl.load_config(argv)
        self.assertEqual(config.a, os.path.expanduser('~/1/2/3'))
        self.assertEqual(config.b, os.path.expanduser('~'))
        self.assertEqual(config.c, os.path.expanduser('~/'))
        self.assertEqual(config.d, '~/')
    
    def test_extra_args(self):
        cl = self.klass(log=log)
        config = cl.load_config(['--a=5', 'b', '--c=10', 'd'])
        self.assertEqual(cl.extra_args, ['b', 'd'])
        self.assertEqual(config.a, 5)
        self.assertEqual(config.c, 10)
        config = cl.load_config(['--', '--a=5', '--c=10'])
        self.assertEqual(cl.extra_args, ['--a=5', '--c=10'])
    
    def test_unicode_args(self):
        cl = self.klass(log=log)
        argv = [u'--a=épsîlön']
        config = cl.load_config(argv)
        self.assertEqual(config.a, u'épsîlön')
    
    def test_unicode_bytes_args(self):
        uarg = u'--a=é'
        try:
            barg = uarg.encode(sys.stdin.encoding)
        except (TypeError, UnicodeEncodeError):
            raise skip("sys.stdin.encoding can't handle 'é'")
        
        cl = self.klass(log=log)
        config = cl.load_config([barg])
        self.assertEqual(config.a, u'é')
    
    def test_unicode_alias(self):
        cl = self.klass(log=log)
        argv = [u'--a=épsîlön']
        config = cl.load_config(argv, aliases=dict(a='A.a'))
        self.assertEqual(config.A.a, u'épsîlön')


class TestArgParseKVCL(TestKeyValueCL):
    klass = KVArgParseConfigLoader

    def test_expanduser2(self):
        cl = self.klass(log=log)
        argv = ['-a', '~/1/2/3', '--b', "'~/1/2/3'"]
        config = cl.load_config(argv, aliases=dict(a='A.a', b='A.b'))
        self.assertEqual(config.A.a, os.path.expanduser('~/1/2/3'))
        self.assertEqual(config.A.b, '~/1/2/3')
    
    def test_eval(self):
        cl = self.klass(log=log)
        argv = ['-c', 'a=5']
        config = cl.load_config(argv, aliases=dict(c='A.c'))
        self.assertEqual(config.A.c, u"a=5")
    

class TestConfig(TestCase):

    def test_setget(self):
        c = Config()
        c.a = 10
        self.assertEqual(c.a, 10)
        self.assertEqual('b' in c, False)

    def test_auto_section(self):
        c = Config()
        self.assertNotIn('A', c)
        assert not c._has_section('A')
        A = c.A
        A.foo = 'hi there'
        self.assertIn('A', c)
        assert c._has_section('A')
        self.assertEqual(c.A.foo, 'hi there')
        del c.A
        self.assertEqual(c.A, Config())

    def test_merge_doesnt_exist(self):
        c1 = Config()
        c2 = Config()
        c2.bar = 10
        c2.Foo.bar = 10
        c1.merge(c2)
        self.assertEqual(c1.Foo.bar, 10)
        self.assertEqual(c1.bar, 10)
        c2.Bar.bar = 10
        c1.merge(c2)
        self.assertEqual(c1.Bar.bar, 10)

    def test_merge_exists(self):
        c1 = Config()
        c2 = Config()
        c1.Foo.bar = 10
        c1.Foo.bam = 30
        c2.Foo.bar = 20
        c2.Foo.wow = 40
        c1.merge(c2)
        self.assertEqual(c1.Foo.bam, 30)
        self.assertEqual(c1.Foo.bar, 20)
        self.assertEqual(c1.Foo.wow, 40)
        c2.Foo.Bam.bam = 10
        c1.merge(c2)
        self.assertEqual(c1.Foo.Bam.bam, 10)

    def test_deepcopy(self):
        c1 = Config()
        c1.Foo.bar = 10
        c1.Foo.bam = 30
        c1.a = 'asdf'
        c1.b = range(10)
        c1.Test.logger = logging.Logger('test')
        c1.Test.get_logger = logging.getLogger('test')
        c2 = copy.deepcopy(c1)
        self.assertEqual(c1, c2)
        self.assertTrue(c1 is not c2)
        self.assertTrue(c1.Foo is not c2.Foo)
        self.assertTrue(c1.Test is not c2.Test)
        self.assertTrue(c1.Test.logger is c2.Test.logger)
        self.assertTrue(c1.Test.get_logger is c2.Test.get_logger)

    def test_builtin(self):
        c1 = Config()
        c1.format = "json"
    
    def test_fromdict(self):
        c1 = Config({'Foo' : {'bar' : 1}})
        self.assertEqual(c1.Foo.__class__, Config)
        self.assertEqual(c1.Foo.bar, 1)
    
    def test_fromdictmerge(self):
        c1 = Config()
        c2 = Config({'Foo' : {'bar' : 1}})
        c1.merge(c2)
        self.assertEqual(c1.Foo.__class__, Config)
        self.assertEqual(c1.Foo.bar, 1)

    def test_fromdictmerge2(self):
        c1 = Config({'Foo' : {'baz' : 2}})
        c2 = Config({'Foo' : {'bar' : 1}})
        c1.merge(c2)
        self.assertEqual(c1.Foo.__class__, Config)
        self.assertEqual(c1.Foo.bar, 1)
        self.assertEqual(c1.Foo.baz, 2)
        self.assertNotIn('baz', c2.Foo)
    
    def test_contains(self):
        c1 = Config({'Foo' : {'baz' : 2}})
        c2 = Config({'Foo' : {'bar' : 1}})
        self.assertIn('Foo', c1)
        self.assertIn('Foo.baz', c1)
        self.assertIn('Foo.bar', c2)
        self.assertNotIn('Foo.bar', c1)
    
    def test_pickle_config(self):
        cfg = Config()
        cfg.Foo.bar = 1
        pcfg = pickle.dumps(cfg)
        cfg2 = pickle.loads(pcfg)
        self.assertEqual(cfg2, cfg)
    
    def test_getattr_section(self):
        cfg = Config()
        self.assertNotIn('Foo', cfg)
        Foo = cfg.Foo
        assert isinstance(Foo, Config)
        self.assertIn('Foo', cfg)

    def test_getitem_section(self):
        cfg = Config()
        self.assertNotIn('Foo', cfg)
        Foo = cfg['Foo']
        assert isinstance(Foo, Config)
        self.assertIn('Foo', cfg)

    def test_getattr_not_section(self):
        cfg = Config()
        self.assertNotIn('foo', cfg)
        foo = cfg.foo
        assert isinstance(foo, LazyConfigValue)
        self.assertIn('foo', cfg)

    def test_getattr_private_missing(self):
        cfg = Config()
        self.assertNotIn('_repr_html_', cfg)
        with self.assertRaises(AttributeError):
            _ = cfg._repr_html_
        self.assertNotIn('_repr_html_', cfg)
        self.assertEqual(len(cfg), 0)

    def test_getitem_not_section(self):
        cfg = Config()
        self.assertNotIn('foo', cfg)
        foo = cfg['foo']
        assert isinstance(foo, LazyConfigValue)
        self.assertIn('foo', cfg)
    
    def test_merge_no_copies(self):
        c = Config()
        c2 = Config()
        c2.Foo.trait = []
        c.merge(c2)
        c2.Foo.trait.append(1)
        self.assertIs(c.Foo, c2.Foo)
        self.assertEqual(c.Foo.trait, [1])
        self.assertEqual(c2.Foo.trait, [1])

