# coding: utf-8
"""
Tests for traitlets.config.application.Application
"""

# Copyright (c) IPython Development Team.
# Distributed under the terms of the Modified BSD License.

import json
import logging
import os
from io import StringIO
from unittest import TestCase

try:
    from unittest import mock
except ImportError:
    import mock

pjoin = os.path.join

from pytest import mark

from traitlets.config.configurable import Configurable
from traitlets.config.loader import Config
from traitlets.tests.utils import check_help_output, check_help_all_output

from traitlets.config.application import (
    Application
)

from ipython_genutils.tempdir import TemporaryDirectory
from traitlets.traitlets import (
    Bool, Unicode, Integer, List, Dict
)


class Foo(Configurable):

    i = Integer(0, help="The integer i.").tag(config=True)
    j = Integer(1, help="The integer j.").tag(config=True)
    name = Unicode(u'Brian', help="First name.").tag(config=True)


class Bar(Configurable):

    b = Integer(0, help="The integer b.").tag(config=True)
    enabled = Bool(True, help="Enable bar.").tag(config=True)


class MyApp(Application):

    name = Unicode(u'myapp')
    running = Bool(False, help="Is the app running?").tag(config=True)
    classes = List([Bar, Foo])
    config_file = Unicode(u'', help="Load this config file").tag(config=True)

    warn_tpyo = Unicode(u"yes the name is wrong on purpose", config=True,
            help="Should print a warning if `MyApp.warn-typo=...` command is passed")

    aliases = Dict({
                    'i' : 'Foo.i',
                    'j' : 'Foo.j',
                    'name' : 'Foo.name',
                    'enabled' : 'Bar.enabled',
                    'log-level' : 'Application.log_level',
                })

    flags = Dict(dict(enable=({'Bar': {'enabled' : True}}, "Set Bar.enabled to True"),
                  disable=({'Bar': {'enabled' : False}}, "Set Bar.enabled to False"),
                  crit=({'Application' : {'log_level' : logging.CRITICAL}},
                        "set level=CRITICAL"),
            ))

    def init_foo(self):
        self.foo = Foo(parent=self)

    def init_bar(self):
        self.bar = Bar(parent=self)


class TestApplication(TestCase):

    def test_log(self):
        stream = StringIO()
        app = MyApp(log_level=logging.INFO)
        handler = logging.StreamHandler(stream)
        # trigger reconstruction of the log formatter
        app.log.handlers = [handler]
        app.log_format = "%(message)s"
        app.log_datefmt = "%Y-%m-%d %H:%M"
        app.log.info("hello")
        assert "hello" in stream.getvalue()

    def test_basic(self):
        app = MyApp()
        self.assertEqual(app.name, u'myapp')
        self.assertEqual(app.running, False)
        self.assertEqual(app.classes, [MyApp,Bar,Foo])
        self.assertEqual(app.config_file, u'')

    def test_config(self):
        app = MyApp()
        app.parse_command_line(["--i=10","--Foo.j=10","--enabled=False","--log-level=50"])
        config = app.config
        self.assertEqual(config.Foo.i, 10)
        self.assertEqual(config.Foo.j, 10)
        self.assertEqual(config.Bar.enabled, False)
        self.assertEqual(config.MyApp.log_level,50)

    def test_config_propagation(self):
        app = MyApp()
        app.parse_command_line(["--i=10","--Foo.j=10","--enabled=False","--log-level=50"])
        app.init_foo()
        app.init_bar()
        self.assertEqual(app.foo.i, 10)
        self.assertEqual(app.foo.j, 10)
        self.assertEqual(app.bar.enabled, False)

    def test_cli_priority(self):
        """Test that loading config files does not override CLI options"""
        name = 'config.py'
        class TestApp(Application):
            value = Unicode().tag(config=True)
            config_file_loaded = Bool().tag(config=True)
            aliases = {'v': 'TestApp.value'}
        app = TestApp()
        with TemporaryDirectory() as td:
            config_file = pjoin(td, name)
            with open(config_file, 'w') as f:
                f.writelines([
                    "c.TestApp.value = 'config file'\n",
                    "c.TestApp.config_file_loaded = True\n"
                ])

            app.parse_command_line(['--v=cli'])
            assert 'value' in app.config.TestApp
            assert app.config.TestApp.value == 'cli'
            assert app.value == 'cli'

            app.load_config_file(name, path=[td])
            assert app.config_file_loaded
            assert app.config.TestApp.value == 'cli'
            assert app.value == 'cli'

    def test_ipython_cli_priority(self):
        # this test is almost entirely redundant with above,
        # but we can keep it around in case of subtle issues creeping into
        # the exact sequence IPython follows.
        name = 'config.py'
        class TestApp(Application):
            value = Unicode().tag(config=True)
            config_file_loaded = Bool().tag(config=True)
            aliases = {'v': 'TestApp.value'}
        app = TestApp()
        with TemporaryDirectory() as td:
            config_file = pjoin(td, name)
            with open(config_file, 'w') as f:
                f.writelines([
                    "c.TestApp.value = 'config file'\n",
                    "c.TestApp.config_file_loaded = True\n"
                ])
            # follow IPython's config-loading sequence to ensure CLI priority is preserved
            app.parse_command_line(['--v=cli'])
            # this is where IPython makes a mistake:
            # it assumes app.config will not be modified,
            # and storing a reference is storing a copy
            cli_config = app.config
            assert 'value' in app.config.TestApp
            assert app.config.TestApp.value == 'cli'
            assert app.value == 'cli'
            app.load_config_file(name, path=[td])
            assert app.config_file_loaded
            # enforce cl-opts override config file opts:
            # this is where IPython makes a mistake: it assumes
            # that cl_config is a different object, but it isn't.
            app.update_config(cli_config)
            assert app.config.TestApp.value == 'cli'
            assert app.value == 'cli'

    def test_flags(self):
        app = MyApp()
        app.parse_command_line(["--disable"])
        app.init_bar()
        self.assertEqual(app.bar.enabled, False)
        app.parse_command_line(["--enable"])
        app.init_bar()
        self.assertEqual(app.bar.enabled, True)

    def test_aliases(self):
        app = MyApp()
        app.parse_command_line(["--i=5", "--j=10"])
        app.init_foo()
        self.assertEqual(app.foo.i, 5)
        app.init_foo()
        self.assertEqual(app.foo.j, 10)

    def test_flag_clobber(self):
        """test that setting flags doesn't clobber existing settings"""
        app = MyApp()
        app.parse_command_line(["--Bar.b=5", "--disable"])
        app.init_bar()
        self.assertEqual(app.bar.enabled, False)
        self.assertEqual(app.bar.b, 5)
        app.parse_command_line(["--enable", "--Bar.b=10"])
        app.init_bar()
        self.assertEqual(app.bar.enabled, True)
        self.assertEqual(app.bar.b, 10)

    def test_warn_autocorrect(self):
        stream = StringIO()
        app = MyApp(log_level=logging.INFO)
        app.log.handlers = [logging.StreamHandler(stream)]

        cfg = Config()
        cfg.MyApp.warn_typo = "WOOOO"
        app.config = cfg

        self.assertIn("warn_typo", stream.getvalue())
        self.assertIn("warn_tpyo", stream.getvalue())


    def test_flatten_flags(self):
        cfg = Config()
        cfg.MyApp.log_level = logging.WARN
        app = MyApp()
        app.update_config(cfg)
        self.assertEqual(app.log_level, logging.WARN)
        self.assertEqual(app.config.MyApp.log_level, logging.WARN)
        app.initialize(["--crit"])
        self.assertEqual(app.log_level, logging.CRITICAL)
        # this would be app.config.Application.log_level if it failed:
        self.assertEqual(app.config.MyApp.log_level, logging.CRITICAL)

    def test_flatten_aliases(self):
        cfg = Config()
        cfg.MyApp.log_level = logging.WARN
        app = MyApp()
        app.update_config(cfg)
        self.assertEqual(app.log_level, logging.WARN)
        self.assertEqual(app.config.MyApp.log_level, logging.WARN)
        app.initialize(["--log-level", "CRITICAL"])
        self.assertEqual(app.log_level, logging.CRITICAL)
        # this would be app.config.Application.log_level if it failed:
        self.assertEqual(app.config.MyApp.log_level, "CRITICAL")

    def test_extra_args(self):
        app = MyApp()
        app.parse_command_line(["--Bar.b=5", 'extra', "--disable", 'args'])
        app.init_bar()
        self.assertEqual(app.bar.enabled, False)
        self.assertEqual(app.bar.b, 5)
        self.assertEqual(app.extra_args, ['extra', 'args'])
        app = MyApp()
        app.parse_command_line(["--Bar.b=5", '--', 'extra', "--disable", 'args'])
        app.init_bar()
        self.assertEqual(app.bar.enabled, True)
        self.assertEqual(app.bar.b, 5)
        self.assertEqual(app.extra_args, ['extra', '--disable', 'args'])

    def test_unicode_argv(self):
        app = MyApp()
        app.parse_command_line(['ünîcødé'])

    def test_document_config_option(self):
        app = MyApp()
        app.document_config_options()

    def test_generate_config_file(self):
        app = MyApp()
        assert 'The integer b.' in app.generate_config_file()

    def test_generate_config_file_classes_to_include(self):
        class NoTraits(Foo, Bar):
            pass

        app = MyApp()
        app.classes.append(NoTraits)
        conf_txt = app.generate_config_file()
        self.assertIn('The integer b.', conf_txt)
        self.assertIn('# Bar(Configurable)', conf_txt)
        self.assertIn('# Foo(Configurable)', conf_txt)
        self.assertNotIn('# Configurable', conf_txt)
        self.assertIn('# NoTraits(Foo,Bar)', conf_txt)

    def test_multi_file(self):
        app = MyApp()
        app.log = logging.getLogger()
        name = 'config.py'
        with TemporaryDirectory('_1') as td1:
            with open(pjoin(td1, name), 'w') as f1:
                f1.write("get_config().MyApp.Bar.b = 1")
            with TemporaryDirectory('_2') as td2:
                with open(pjoin(td2, name), 'w') as f2:
                    f2.write("get_config().MyApp.Bar.b = 2")
                app.load_config_file(name, path=[td2, td1])
                app.init_bar()
                self.assertEqual(app.bar.b, 2)
                app.load_config_file(name, path=[td1, td2])
                app.init_bar()
                self.assertEqual(app.bar.b, 1)

    @mark.skipif(not hasattr(TestCase, 'assertLogs'), reason='requires TestCase.assertLogs')
    def test_log_collisions(self):
        app = MyApp()
        app.log = logging.getLogger()
        app.log.setLevel(logging.INFO)
        name = 'config'
        with TemporaryDirectory('_1') as td:
            with open(pjoin(td, name + '.py'), 'w') as f:
                f.write("get_config().Bar.b = 1")
            with open(pjoin(td, name + '.json'), 'w') as f:
                json.dump({
                    'Bar': {
                        'b': 2
                    }
                }, f)
            with self.assertLogs(app.log, logging.WARNING) as captured:
                app.load_config_file(name, path=[td])
                app.init_bar()
        assert app.bar.b == 2
        output = '\n'.join(captured.output)
        assert 'Collision' in output
        assert '1 ignored, using 2' in output
        assert pjoin(td, name + '.py') in output
        assert pjoin(td, name + '.json') in output

    @mark.skipif(not hasattr(TestCase, 'assertLogs'), reason='requires TestCase.assertLogs')
    def test_log_bad_config(self):
        app = MyApp()
        app.log = logging.getLogger()
        name = 'config.py'
        with TemporaryDirectory() as td:
            with open(pjoin(td, name), 'w') as f:
                f.write("syntax error()")
            with self.assertLogs(app.log, logging.ERROR) as captured:
                app.load_config_file(name, path=[td])
        output = '\n'.join(captured.output)
        self.assertIn('SyntaxError', output)

    def test_raise_on_bad_config(self):
        app = MyApp()
        app.raise_config_file_errors = True
        app.log = logging.getLogger()
        name = 'config.py'
        with TemporaryDirectory() as td:
            with open(pjoin(td, name), 'w') as f:
                f.write("syntax error()")
            with self.assertRaises(SyntaxError):
                app.load_config_file(name, path=[td])


class DeprecatedApp(Application):
    override_called = False
    parent_called = False
    def _config_changed(self, name, old, new):
        self.override_called = True
        def _capture(*args):
            self.parent_called = True
        with mock.patch.object(self.log, 'debug', _capture):
            super(DeprecatedApp, self)._config_changed(name, old, new)


def test_deprecated_notifier():
    app = DeprecatedApp()
    assert not app.override_called
    assert not app.parent_called
    app.config = Config({'A': {'b': 'c'}})
    assert app.override_called
    assert app.parent_called


def test_help_output():
    check_help_output(__name__)
    check_help_all_output(__name__)

if __name__ == '__main__':
    # for test_help_output:
    MyApp.launch_instance()