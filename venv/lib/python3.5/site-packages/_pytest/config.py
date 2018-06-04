""" command line options, ini-file and conftest.py processing. """
from __future__ import absolute_import, division, print_function
import argparse
import shlex
import traceback
import types
import warnings
import copy
import six
import py
# DON't import pytest here because it causes import cycle troubles
import sys
import os
from _pytest.outcomes import Skipped

import _pytest._code
import _pytest.hookspec  # the extension point definitions
import _pytest.assertion
from pluggy import PluginManager, HookimplMarker, HookspecMarker
from _pytest.compat import safe_str

hookimpl = HookimplMarker("pytest")
hookspec = HookspecMarker("pytest")

# pytest startup
#


class ConftestImportFailure(Exception):
    def __init__(self, path, excinfo):
        Exception.__init__(self, path, excinfo)
        self.path = path
        self.excinfo = excinfo

    def __str__(self):
        etype, evalue, etb = self.excinfo
        formatted = traceback.format_tb(etb)
        # The level of the tracebacks we want to print is hand crafted :(
        return repr(evalue) + '\n' + ''.join(formatted[2:])


def main(args=None, plugins=None):
    """ return exit code, after performing an in-process test run.

    :arg args: list of command line arguments.

    :arg plugins: list of plugin objects to be auto-registered during
                  initialization.
    """
    try:
        try:
            config = _prepareconfig(args, plugins)
        except ConftestImportFailure as e:
            tw = py.io.TerminalWriter(sys.stderr)
            for line in traceback.format_exception(*e.excinfo):
                tw.line(line.rstrip(), red=True)
            tw.line("ERROR: could not load %s\n" % (e.path,), red=True)
            return 4
        else:
            try:
                return config.hook.pytest_cmdline_main(config=config)
            finally:
                config._ensure_unconfigure()
    except UsageError as e:
        tw = py.io.TerminalWriter(sys.stderr)
        for msg in e.args:
            tw.line("ERROR: {}\n".format(msg), red=True)
        return 4


class cmdline(object):  # NOQA compatibility namespace
    main = staticmethod(main)


class UsageError(Exception):
    """ error in pytest usage or invocation"""


class PrintHelp(Exception):
    """Raised when pytest should print it's help to skip the rest of the
    argument parsing and validation."""
    pass


def filename_arg(path, optname):
    """ Argparse type validator for filename arguments.

    :path: path of filename
    :optname: name of the option
    """
    if os.path.isdir(path):
        raise UsageError("{0} must be a filename, given: {1}".format(optname, path))
    return path


def directory_arg(path, optname):
    """Argparse type validator for directory arguments.

    :path: path of directory
    :optname: name of the option
    """
    if not os.path.isdir(path):
        raise UsageError("{0} must be a directory, given: {1}".format(optname, path))
    return path


default_plugins = (
    "mark main terminal runner python fixtures debugging unittest capture skipping "
    "tmpdir monkeypatch recwarn pastebin helpconfig nose assertion "
    "junitxml resultlog doctest cacheprovider freeze_support "
    "setuponly setupplan warnings logging").split()


builtin_plugins = set(default_plugins)
builtin_plugins.add("pytester")


def get_config():
    # subsequent calls to main will create a fresh instance
    pluginmanager = PytestPluginManager()
    config = Config(pluginmanager)
    for spec in default_plugins:
        pluginmanager.import_plugin(spec)
    return config


def get_plugin_manager():
    """
    Obtain a new instance of the
    :py:class:`_pytest.config.PytestPluginManager`, with default plugins
    already loaded.

    This function can be used by integration with other tools, like hooking
    into pytest to run tests into an IDE.
    """
    return get_config().pluginmanager


def _prepareconfig(args=None, plugins=None):
    warning = None
    if args is None:
        args = sys.argv[1:]
    elif isinstance(args, py.path.local):
        args = [str(args)]
    elif not isinstance(args, (tuple, list)):
        if not isinstance(args, str):
            raise ValueError("not a string or argument list: %r" % (args,))
        args = shlex.split(args, posix=sys.platform != "win32")
        from _pytest import deprecated
        warning = deprecated.MAIN_STR_ARGS
    config = get_config()
    pluginmanager = config.pluginmanager
    try:
        if plugins:
            for plugin in plugins:
                if isinstance(plugin, six.string_types):
                    pluginmanager.consider_pluginarg(plugin)
                else:
                    pluginmanager.register(plugin)
        if warning:
            config.warn('C1', warning)
        return pluginmanager.hook.pytest_cmdline_parse(
            pluginmanager=pluginmanager, args=args)
    except BaseException:
        config._ensure_unconfigure()
        raise


class PytestPluginManager(PluginManager):
    """
    Overwrites :py:class:`pluggy.PluginManager <pluggy.PluginManager>` to add pytest-specific
    functionality:

    * loading plugins from the command line, ``PYTEST_PLUGINS`` env variable and
      ``pytest_plugins`` global variables found in plugins being loaded;
    * ``conftest.py`` loading during start-up;
    """

    def __init__(self):
        super(PytestPluginManager, self).__init__("pytest", implprefix="pytest_")
        self._conftest_plugins = set()

        # state related to local conftest plugins
        self._path2confmods = {}
        self._conftestpath2mod = {}
        self._confcutdir = None
        self._noconftest = False
        self._duplicatepaths = set()

        self.add_hookspecs(_pytest.hookspec)
        self.register(self)
        if os.environ.get('PYTEST_DEBUG'):
            err = sys.stderr
            encoding = getattr(err, 'encoding', 'utf8')
            try:
                err = py.io.dupfile(err, encoding=encoding)
            except Exception:
                pass
            self.trace.root.setwriter(err.write)
            self.enable_tracing()

        # Config._consider_importhook will set a real object if required.
        self.rewrite_hook = _pytest.assertion.DummyRewriteHook()
        # Used to know when we are importing conftests after the pytest_configure stage
        self._configured = False

    def addhooks(self, module_or_class):
        """
        .. deprecated:: 2.8

        Use :py:meth:`pluggy.PluginManager.add_hookspecs <PluginManager.add_hookspecs>`
        instead.
        """
        warning = dict(code="I2",
                       fslocation=_pytest._code.getfslineno(sys._getframe(1)),
                       nodeid=None,
                       message="use pluginmanager.add_hookspecs instead of "
                               "deprecated addhooks() method.")
        self._warn(warning)
        return self.add_hookspecs(module_or_class)

    def parse_hookimpl_opts(self, plugin, name):
        # pytest hooks are always prefixed with pytest_
        # so we avoid accessing possibly non-readable attributes
        # (see issue #1073)
        if not name.startswith("pytest_"):
            return
        # ignore some historic special names which can not be hooks anyway
        if name == "pytest_plugins" or name.startswith("pytest_funcarg__"):
            return

        method = getattr(plugin, name)
        opts = super(PytestPluginManager, self).parse_hookimpl_opts(plugin, name)
        if opts is not None:
            for name in ("tryfirst", "trylast", "optionalhook", "hookwrapper"):
                opts.setdefault(name, hasattr(method, name))
        return opts

    def parse_hookspec_opts(self, module_or_class, name):
        opts = super(PytestPluginManager, self).parse_hookspec_opts(
            module_or_class, name)
        if opts is None:
            method = getattr(module_or_class, name)
            if name.startswith("pytest_"):
                opts = {"firstresult": hasattr(method, "firstresult"),
                        "historic": hasattr(method, "historic")}
        return opts

    def register(self, plugin, name=None):
        if name in ['pytest_catchlog', 'pytest_capturelog']:
            self._warn('{0} plugin has been merged into the core, '
                       'please remove it from your requirements.'.format(
                           name.replace('_', '-')))
            return
        ret = super(PytestPluginManager, self).register(plugin, name)
        if ret:
            self.hook.pytest_plugin_registered.call_historic(
                kwargs=dict(plugin=plugin, manager=self))

            if isinstance(plugin, types.ModuleType):
                self.consider_module(plugin)
        return ret

    def getplugin(self, name):
        # support deprecated naming because plugins (xdist e.g.) use it
        return self.get_plugin(name)

    def hasplugin(self, name):
        """Return True if the plugin with the given name is registered."""
        return bool(self.get_plugin(name))

    def pytest_configure(self, config):
        # XXX now that the pluginmanager exposes hookimpl(tryfirst...)
        # we should remove tryfirst/trylast as markers
        config.addinivalue_line("markers",
                                "tryfirst: mark a hook implementation function such that the "
                                "plugin machinery will try to call it first/as early as possible.")
        config.addinivalue_line("markers",
                                "trylast: mark a hook implementation function such that the "
                                "plugin machinery will try to call it last/as late as possible.")
        self._configured = True

    def _warn(self, message):
        kwargs = message if isinstance(message, dict) else {
            'code': 'I1',
            'message': message,
            'fslocation': None,
            'nodeid': None,
        }
        self.hook.pytest_logwarning.call_historic(kwargs=kwargs)

    #
    # internal API for local conftest plugin handling
    #
    def _set_initial_conftests(self, namespace):
        """ load initial conftest files given a preparsed "namespace".
            As conftest files may add their own command line options
            which have arguments ('--my-opt somepath') we might get some
            false positives.  All builtin and 3rd party plugins will have
            been loaded, however, so common options will not confuse our logic
            here.
        """
        current = py.path.local()
        self._confcutdir = current.join(namespace.confcutdir, abs=True) \
            if namespace.confcutdir else None
        self._noconftest = namespace.noconftest
        testpaths = namespace.file_or_dir
        foundanchor = False
        for path in testpaths:
            path = str(path)
            # remove node-id syntax
            i = path.find("::")
            if i != -1:
                path = path[:i]
            anchor = current.join(path, abs=1)
            if exists(anchor):  # we found some file object
                self._try_load_conftest(anchor)
                foundanchor = True
        if not foundanchor:
            self._try_load_conftest(current)

    def _try_load_conftest(self, anchor):
        self._getconftestmodules(anchor)
        # let's also consider test* subdirs
        if anchor.check(dir=1):
            for x in anchor.listdir("test*"):
                if x.check(dir=1):
                    self._getconftestmodules(x)

    def _getconftestmodules(self, path):
        if self._noconftest:
            return []
        try:
            return self._path2confmods[path]
        except KeyError:
            if path.isfile():
                clist = self._getconftestmodules(path.dirpath())
            else:
                # XXX these days we may rather want to use config.rootdir
                # and allow users to opt into looking into the rootdir parent
                # directories instead of requiring to specify confcutdir
                clist = []
                for parent in path.parts():
                    if self._confcutdir and self._confcutdir.relto(parent):
                        continue
                    conftestpath = parent.join("conftest.py")
                    if conftestpath.isfile():
                        mod = self._importconftest(conftestpath)
                        clist.append(mod)

            self._path2confmods[path] = clist
            return clist

    def _rget_with_confmod(self, name, path):
        modules = self._getconftestmodules(path)
        for mod in reversed(modules):
            try:
                return mod, getattr(mod, name)
            except AttributeError:
                continue
        raise KeyError(name)

    def _importconftest(self, conftestpath):
        try:
            return self._conftestpath2mod[conftestpath]
        except KeyError:
            pkgpath = conftestpath.pypkgpath()
            if pkgpath is None:
                _ensure_removed_sysmodule(conftestpath.purebasename)
            try:
                mod = conftestpath.pyimport()
                if hasattr(mod, 'pytest_plugins') and self._configured:
                    from _pytest.deprecated import PYTEST_PLUGINS_FROM_NON_TOP_LEVEL_CONFTEST
                    warnings.warn(PYTEST_PLUGINS_FROM_NON_TOP_LEVEL_CONFTEST)
            except Exception:
                raise ConftestImportFailure(conftestpath, sys.exc_info())

            self._conftest_plugins.add(mod)
            self._conftestpath2mod[conftestpath] = mod
            dirpath = conftestpath.dirpath()
            if dirpath in self._path2confmods:
                for path, mods in self._path2confmods.items():
                    if path and path.relto(dirpath) or path == dirpath:
                        assert mod not in mods
                        mods.append(mod)
            self.trace("loaded conftestmodule %r" % (mod))
            self.consider_conftest(mod)
            return mod

    #
    # API for bootstrapping plugin loading
    #
    #

    def consider_preparse(self, args):
        for opt1, opt2 in zip(args, args[1:]):
            if opt1 == "-p":
                self.consider_pluginarg(opt2)

    def consider_pluginarg(self, arg):
        if arg.startswith("no:"):
            name = arg[3:]
            self.set_blocked(name)
            if not name.startswith("pytest_"):
                self.set_blocked("pytest_" + name)
        else:
            self.import_plugin(arg)

    def consider_conftest(self, conftestmodule):
        self.register(conftestmodule, name=conftestmodule.__file__)

    def consider_env(self):
        self._import_plugin_specs(os.environ.get("PYTEST_PLUGINS"))

    def consider_module(self, mod):
        self._import_plugin_specs(getattr(mod, 'pytest_plugins', []))

    def _import_plugin_specs(self, spec):
        plugins = _get_plugin_specs_as_list(spec)
        for import_spec in plugins:
            self.import_plugin(import_spec)

    def import_plugin(self, modname):
        # most often modname refers to builtin modules, e.g. "pytester",
        # "terminal" or "capture".  Those plugins are registered under their
        # basename for historic purposes but must be imported with the
        # _pytest prefix.
        assert isinstance(modname, (six.text_type, str)), "module name as text required, got %r" % modname
        modname = str(modname)
        if self.is_blocked(modname) or self.get_plugin(modname) is not None:
            return
        if modname in builtin_plugins:
            importspec = "_pytest." + modname
        else:
            importspec = modname
        self.rewrite_hook.mark_rewrite(importspec)
        try:
            __import__(importspec)
        except ImportError as e:
            new_exc_type = ImportError
            new_exc_message = 'Error importing plugin "%s": %s' % (modname, safe_str(e.args[0]))
            new_exc = new_exc_type(new_exc_message)

            six.reraise(new_exc_type, new_exc, sys.exc_info()[2])

        except Skipped as e:
            self._warn("skipped plugin %r: %s" % ((modname, e.msg)))
        else:
            mod = sys.modules[importspec]
            self.register(mod, modname)


def _get_plugin_specs_as_list(specs):
    """
    Parses a list of "plugin specs" and returns a list of plugin names.

    Plugin specs can be given as a list of strings separated by "," or already as a list/tuple in
    which case it is returned as a list. Specs can also be `None` in which case an
    empty list is returned.
    """
    if specs is not None:
        if isinstance(specs, str):
            specs = specs.split(',') if specs else []
        if not isinstance(specs, (list, tuple)):
            raise UsageError("Plugin specs must be a ','-separated string or a "
                             "list/tuple of strings for plugin names. Given: %r" % specs)
        return list(specs)
    return []


class Parser(object):
    """ Parser for command line arguments and ini-file values.

    :ivar extra_info: dict of generic param -> value to display in case
        there's an error processing the command line arguments.
    """

    def __init__(self, usage=None, processopt=None):
        self._anonymous = OptionGroup("custom options", parser=self)
        self._groups = []
        self._processopt = processopt
        self._usage = usage
        self._inidict = {}
        self._ininames = []
        self.extra_info = {}

    def processoption(self, option):
        if self._processopt:
            if option.dest:
                self._processopt(option)

    def getgroup(self, name, description="", after=None):
        """ get (or create) a named option Group.

        :name: name of the option group.
        :description: long description for --help output.
        :after: name of other group, used for ordering --help output.

        The returned group object has an ``addoption`` method with the same
        signature as :py:func:`parser.addoption
        <_pytest.config.Parser.addoption>` but will be shown in the
        respective group in the output of ``pytest. --help``.
        """
        for group in self._groups:
            if group.name == name:
                return group
        group = OptionGroup(name, description, parser=self)
        i = 0
        for i, grp in enumerate(self._groups):
            if grp.name == after:
                break
        self._groups.insert(i + 1, group)
        return group

    def addoption(self, *opts, **attrs):
        """ register a command line option.

        :opts: option names, can be short or long options.
        :attrs: same attributes which the ``add_option()`` function of the
           `argparse library
           <http://docs.python.org/2/library/argparse.html>`_
           accepts.

        After command line parsing options are available on the pytest config
        object via ``config.option.NAME`` where ``NAME`` is usually set
        by passing a ``dest`` attribute, for example
        ``addoption("--long", dest="NAME", ...)``.
        """
        self._anonymous.addoption(*opts, **attrs)

    def parse(self, args, namespace=None):
        from _pytest._argcomplete import try_argcomplete
        self.optparser = self._getparser()
        try_argcomplete(self.optparser)
        return self.optparser.parse_args([str(x) for x in args], namespace=namespace)

    def _getparser(self):
        from _pytest._argcomplete import filescompleter
        optparser = MyOptionParser(self, self.extra_info)
        groups = self._groups + [self._anonymous]
        for group in groups:
            if group.options:
                desc = group.description or group.name
                arggroup = optparser.add_argument_group(desc)
                for option in group.options:
                    n = option.names()
                    a = option.attrs()
                    arggroup.add_argument(*n, **a)
        # bash like autocompletion for dirs (appending '/')
        optparser.add_argument(FILE_OR_DIR, nargs='*').completer = filescompleter
        return optparser

    def parse_setoption(self, args, option, namespace=None):
        parsedoption = self.parse(args, namespace=namespace)
        for name, value in parsedoption.__dict__.items():
            setattr(option, name, value)
        return getattr(parsedoption, FILE_OR_DIR)

    def parse_known_args(self, args, namespace=None):
        """parses and returns a namespace object with known arguments at this
        point.
        """
        return self.parse_known_and_unknown_args(args, namespace=namespace)[0]

    def parse_known_and_unknown_args(self, args, namespace=None):
        """parses and returns a namespace object with known arguments, and
        the remaining arguments unknown at this point.
        """
        optparser = self._getparser()
        args = [str(x) for x in args]
        return optparser.parse_known_args(args, namespace=namespace)

    def addini(self, name, help, type=None, default=None):
        """ register an ini-file option.

        :name: name of the ini-variable
        :type: type of the variable, can be ``pathlist``, ``args``, ``linelist``
               or ``bool``.
        :default: default value if no ini-file option exists but is queried.

        The value of ini-variables can be retrieved via a call to
        :py:func:`config.getini(name) <_pytest.config.Config.getini>`.
        """
        assert type in (None, "pathlist", "args", "linelist", "bool")
        self._inidict[name] = (help, type, default)
        self._ininames.append(name)


class ArgumentError(Exception):
    """
    Raised if an Argument instance is created with invalid or
    inconsistent arguments.
    """

    def __init__(self, msg, option):
        self.msg = msg
        self.option_id = str(option)

    def __str__(self):
        if self.option_id:
            return "option %s: %s" % (self.option_id, self.msg)
        else:
            return self.msg


class Argument(object):
    """class that mimics the necessary behaviour of optparse.Option

    its currently a least effort implementation
    and ignoring choices and integer prefixes
    https://docs.python.org/3/library/optparse.html#optparse-standard-option-types
    """
    _typ_map = {
        'int': int,
        'string': str,
        'float': float,
        'complex': complex,
    }

    def __init__(self, *names, **attrs):
        """store parms in private vars for use in add_argument"""
        self._attrs = attrs
        self._short_opts = []
        self._long_opts = []
        self.dest = attrs.get('dest')
        if '%default' in (attrs.get('help') or ''):
            warnings.warn(
                'pytest now uses argparse. "%default" should be'
                ' changed to "%(default)s" ',
                DeprecationWarning,
                stacklevel=3)
        try:
            typ = attrs['type']
        except KeyError:
            pass
        else:
            # this might raise a keyerror as well, don't want to catch that
            if isinstance(typ, six.string_types):
                if typ == 'choice':
                    warnings.warn(
                        'type argument to addoption() is a string %r.'
                        ' For parsearg this is optional and when supplied'
                        ' should be a type.'
                        ' (options: %s)' % (typ, names),
                        DeprecationWarning,
                        stacklevel=3)
                    # argparse expects a type here take it from
                    # the type of the first element
                    attrs['type'] = type(attrs['choices'][0])
                else:
                    warnings.warn(
                        'type argument to addoption() is a string %r.'
                        ' For parsearg this should be a type.'
                        ' (options: %s)' % (typ, names),
                        DeprecationWarning,
                        stacklevel=3)
                    attrs['type'] = Argument._typ_map[typ]
                # used in test_parseopt -> test_parse_defaultgetter
                self.type = attrs['type']
            else:
                self.type = typ
        try:
            # attribute existence is tested in Config._processopt
            self.default = attrs['default']
        except KeyError:
            pass
        self._set_opt_strings(names)
        if not self.dest:
            if self._long_opts:
                self.dest = self._long_opts[0][2:].replace('-', '_')
            else:
                try:
                    self.dest = self._short_opts[0][1:]
                except IndexError:
                    raise ArgumentError(
                        'need a long or short option', self)

    def names(self):
        return self._short_opts + self._long_opts

    def attrs(self):
        # update any attributes set by processopt
        attrs = 'default dest help'.split()
        if self.dest:
            attrs.append(self.dest)
        for attr in attrs:
            try:
                self._attrs[attr] = getattr(self, attr)
            except AttributeError:
                pass
        if self._attrs.get('help'):
            a = self._attrs['help']
            a = a.replace('%default', '%(default)s')
            # a = a.replace('%prog', '%(prog)s')
            self._attrs['help'] = a
        return self._attrs

    def _set_opt_strings(self, opts):
        """directly from optparse

        might not be necessary as this is passed to argparse later on"""
        for opt in opts:
            if len(opt) < 2:
                raise ArgumentError(
                    "invalid option string %r: "
                    "must be at least two characters long" % opt, self)
            elif len(opt) == 2:
                if not (opt[0] == "-" and opt[1] != "-"):
                    raise ArgumentError(
                        "invalid short option string %r: "
                        "must be of the form -x, (x any non-dash char)" % opt,
                        self)
                self._short_opts.append(opt)
            else:
                if not (opt[0:2] == "--" and opt[2] != "-"):
                    raise ArgumentError(
                        "invalid long option string %r: "
                        "must start with --, followed by non-dash" % opt,
                        self)
                self._long_opts.append(opt)

    def __repr__(self):
        args = []
        if self._short_opts:
            args += ['_short_opts: ' + repr(self._short_opts)]
        if self._long_opts:
            args += ['_long_opts: ' + repr(self._long_opts)]
        args += ['dest: ' + repr(self.dest)]
        if hasattr(self, 'type'):
            args += ['type: ' + repr(self.type)]
        if hasattr(self, 'default'):
            args += ['default: ' + repr(self.default)]
        return 'Argument({0})'.format(', '.join(args))


class OptionGroup(object):
    def __init__(self, name, description="", parser=None):
        self.name = name
        self.description = description
        self.options = []
        self.parser = parser

    def addoption(self, *optnames, **attrs):
        """ add an option to this group.

        if a shortened version of a long option is specified it will
        be suppressed in the help. addoption('--twowords', '--two-words')
        results in help showing '--two-words' only, but --twowords gets
        accepted **and** the automatic destination is in args.twowords
        """
        conflict = set(optnames).intersection(
            name for opt in self.options for name in opt.names())
        if conflict:
            raise ValueError("option names %s already added" % conflict)
        option = Argument(*optnames, **attrs)
        self._addoption_instance(option, shortupper=False)

    def _addoption(self, *optnames, **attrs):
        option = Argument(*optnames, **attrs)
        self._addoption_instance(option, shortupper=True)

    def _addoption_instance(self, option, shortupper=False):
        if not shortupper:
            for opt in option._short_opts:
                if opt[0] == '-' and opt[1].islower():
                    raise ValueError("lowercase shortoptions reserved")
        if self.parser:
            self.parser.processoption(option)
        self.options.append(option)


class MyOptionParser(argparse.ArgumentParser):
    def __init__(self, parser, extra_info=None):
        if not extra_info:
            extra_info = {}
        self._parser = parser
        argparse.ArgumentParser.__init__(self, usage=parser._usage,
                                         add_help=False, formatter_class=DropShorterLongHelpFormatter)
        # extra_info is a dict of (param -> value) to display if there's
        # an usage error to provide more contextual information to the user
        self.extra_info = extra_info

    def parse_args(self, args=None, namespace=None):
        """allow splitting of positional arguments"""
        args, argv = self.parse_known_args(args, namespace)
        if argv:
            for arg in argv:
                if arg and arg[0] == '-':
                    lines = ['unrecognized arguments: %s' % (' '.join(argv))]
                    for k, v in sorted(self.extra_info.items()):
                        lines.append('  %s: %s' % (k, v))
                    self.error('\n'.join(lines))
            getattr(args, FILE_OR_DIR).extend(argv)
        return args


class DropShorterLongHelpFormatter(argparse.HelpFormatter):
    """shorten help for long options that differ only in extra hyphens

    - collapse **long** options that are the same except for extra hyphens
    - special action attribute map_long_option allows surpressing additional
      long options
    - shortcut if there are only two options and one of them is a short one
    - cache result on action object as this is called at least 2 times
    """

    def _format_action_invocation(self, action):
        orgstr = argparse.HelpFormatter._format_action_invocation(self, action)
        if orgstr and orgstr[0] != '-':  # only optional arguments
            return orgstr
        res = getattr(action, '_formatted_action_invocation', None)
        if res:
            return res
        options = orgstr.split(', ')
        if len(options) == 2 and (len(options[0]) == 2 or len(options[1]) == 2):
            # a shortcut for '-h, --help' or '--abc', '-a'
            action._formatted_action_invocation = orgstr
            return orgstr
        return_list = []
        option_map = getattr(action, 'map_long_option', {})
        if option_map is None:
            option_map = {}
        short_long = {}
        for option in options:
            if len(option) == 2 or option[2] == ' ':
                continue
            if not option.startswith('--'):
                raise ArgumentError('long optional argument without "--": [%s]'
                                    % (option), self)
            xxoption = option[2:]
            if xxoption.split()[0] not in option_map:
                shortened = xxoption.replace('-', '')
                if shortened not in short_long or \
                   len(short_long[shortened]) < len(xxoption):
                    short_long[shortened] = xxoption
        # now short_long has been filled out to the longest with dashes
        # **and** we keep the right option ordering from add_argument
        for option in options:
            if len(option) == 2 or option[2] == ' ':
                return_list.append(option)
            if option[2:] == short_long.get(option.replace('-', '')):
                return_list.append(option.replace(' ', '=', 1))
        action._formatted_action_invocation = ', '.join(return_list)
        return action._formatted_action_invocation


def _ensure_removed_sysmodule(modname):
    try:
        del sys.modules[modname]
    except KeyError:
        pass


class Notset(object):
    def __repr__(self):
        return "<NOTSET>"


notset = Notset()
FILE_OR_DIR = 'file_or_dir'


def _iter_rewritable_modules(package_files):
    for fn in package_files:
        is_simple_module = '/' not in fn and fn.endswith('.py')
        is_package = fn.count('/') == 1 and fn.endswith('__init__.py')
        if is_simple_module:
            module_name, _ = os.path.splitext(fn)
            yield module_name
        elif is_package:
            package_name = os.path.dirname(fn)
            yield package_name


class Config(object):
    """ access to configuration values, pluginmanager and plugin hooks.  """

    def __init__(self, pluginmanager):
        #: access to command line option as attributes.
        #: (deprecated), use :py:func:`getoption() <_pytest.config.Config.getoption>` instead
        self.option = argparse.Namespace()
        _a = FILE_OR_DIR
        self._parser = Parser(
            usage="%%(prog)s [options] [%s] [%s] [...]" % (_a, _a),
            processopt=self._processopt,
        )
        #: a pluginmanager instance
        self.pluginmanager = pluginmanager
        self.trace = self.pluginmanager.trace.root.get("config")
        self.hook = self.pluginmanager.hook
        self._inicache = {}
        self._override_ini = ()
        self._opt2dest = {}
        self._cleanup = []
        self._warn = self.pluginmanager._warn
        self.pluginmanager.register(self, "pytestconfig")
        self._configured = False

        def do_setns(dic):
            import pytest
            setns(pytest, dic)

        self.hook.pytest_namespace.call_historic(do_setns, {})
        self.hook.pytest_addoption.call_historic(kwargs=dict(parser=self._parser))

    def add_cleanup(self, func):
        """ Add a function to be called when the config object gets out of
        use (usually coninciding with pytest_unconfigure)."""
        self._cleanup.append(func)

    def _do_configure(self):
        assert not self._configured
        self._configured = True
        self.hook.pytest_configure.call_historic(kwargs=dict(config=self))

    def _ensure_unconfigure(self):
        if self._configured:
            self._configured = False
            self.hook.pytest_unconfigure(config=self)
            self.hook.pytest_configure._call_history = []
        while self._cleanup:
            fin = self._cleanup.pop()
            fin()

    def warn(self, code, message, fslocation=None, nodeid=None):
        """ generate a warning for this test session. """
        self.hook.pytest_logwarning.call_historic(kwargs=dict(
            code=code, message=message,
            fslocation=fslocation, nodeid=nodeid))

    def get_terminal_writer(self):
        return self.pluginmanager.get_plugin("terminalreporter")._tw

    def pytest_cmdline_parse(self, pluginmanager, args):
        # REF1 assert self == pluginmanager.config, (self, pluginmanager.config)
        self.parse(args)
        return self

    def notify_exception(self, excinfo, option=None):
        if option and option.fulltrace:
            style = "long"
        else:
            style = "native"
        excrepr = excinfo.getrepr(funcargs=True,
                                  showlocals=getattr(option, 'showlocals', False),
                                  style=style,
                                  )
        res = self.hook.pytest_internalerror(excrepr=excrepr,
                                             excinfo=excinfo)
        if not any(res):
            for line in str(excrepr).split("\n"):
                sys.stderr.write("INTERNALERROR> %s\n" % line)
                sys.stderr.flush()

    def cwd_relative_nodeid(self, nodeid):
        # nodeid's are relative to the rootpath, compute relative to cwd
        if self.invocation_dir != self.rootdir:
            fullpath = self.rootdir.join(nodeid)
            nodeid = self.invocation_dir.bestrelpath(fullpath)
        return nodeid

    @classmethod
    def fromdictargs(cls, option_dict, args):
        """ constructor useable for subprocesses. """
        config = get_config()
        config.option.__dict__.update(option_dict)
        config.parse(args, addopts=False)
        for x in config.option.plugins:
            config.pluginmanager.consider_pluginarg(x)
        return config

    def _processopt(self, opt):
        for name in opt._short_opts + opt._long_opts:
            self._opt2dest[name] = opt.dest

        if hasattr(opt, 'default') and opt.dest:
            if not hasattr(self.option, opt.dest):
                setattr(self.option, opt.dest, opt.default)

    @hookimpl(trylast=True)
    def pytest_load_initial_conftests(self, early_config):
        self.pluginmanager._set_initial_conftests(early_config.known_args_namespace)

    def _initini(self, args):
        ns, unknown_args = self._parser.parse_known_and_unknown_args(args, namespace=copy.copy(self.option))
        r = determine_setup(ns.inifilename, ns.file_or_dir + unknown_args, warnfunc=self.warn,
                            rootdir_cmd_arg=ns.rootdir or None)
        self.rootdir, self.inifile, self.inicfg = r
        self._parser.extra_info['rootdir'] = self.rootdir
        self._parser.extra_info['inifile'] = self.inifile
        self.invocation_dir = py.path.local()
        self._parser.addini('addopts', 'extra command line options', 'args')
        self._parser.addini('minversion', 'minimally required pytest version')
        self._override_ini = ns.override_ini or ()

    def _consider_importhook(self, args):
        """Install the PEP 302 import hook if using assertion rewriting.

        Needs to parse the --assert=<mode> option from the commandline
        and find all the installed plugins to mark them for rewriting
        by the importhook.
        """
        ns, unknown_args = self._parser.parse_known_and_unknown_args(args)
        mode = ns.assertmode
        if mode == 'rewrite':
            try:
                hook = _pytest.assertion.install_importhook(self)
            except SystemError:
                mode = 'plain'
            else:
                self._mark_plugins_for_rewrite(hook)
        _warn_about_missing_assertion(mode)

    def _mark_plugins_for_rewrite(self, hook):
        """
        Given an importhook, mark for rewrite any top-level
        modules or packages in the distribution package for
        all pytest plugins.
        """
        import pkg_resources
        self.pluginmanager.rewrite_hook = hook

        # 'RECORD' available for plugins installed normally (pip install)
        # 'SOURCES.txt' available for plugins installed in dev mode (pip install -e)
        # for installed plugins 'SOURCES.txt' returns an empty list, and vice-versa
        # so it shouldn't be an issue
        metadata_files = 'RECORD', 'SOURCES.txt'

        package_files = (
            entry.split(',')[0]
            for entrypoint in pkg_resources.iter_entry_points('pytest11')
            for metadata in metadata_files
            for entry in entrypoint.dist._get_metadata(metadata)
        )

        for name in _iter_rewritable_modules(package_files):
            hook.mark_rewrite(name)

    def _preparse(self, args, addopts=True):
        if addopts:
            args[:] = shlex.split(os.environ.get('PYTEST_ADDOPTS', '')) + args
        self._initini(args)
        if addopts:
            args[:] = self.getini("addopts") + args
        self._checkversion()
        self._consider_importhook(args)
        self.pluginmanager.consider_preparse(args)
        self.pluginmanager.load_setuptools_entrypoints('pytest11')
        self.pluginmanager.consider_env()
        self.known_args_namespace = ns = self._parser.parse_known_args(
            args, namespace=copy.copy(self.option))
        if self.known_args_namespace.confcutdir is None and self.inifile:
            confcutdir = py.path.local(self.inifile).dirname
            self.known_args_namespace.confcutdir = confcutdir
        try:
            self.hook.pytest_load_initial_conftests(early_config=self,
                                                    args=args, parser=self._parser)
        except ConftestImportFailure:
            e = sys.exc_info()[1]
            if ns.help or ns.version:
                # we don't want to prevent --help/--version to work
                # so just let is pass and print a warning at the end
                self._warn("could not load initial conftests (%s)\n" % e.path)
            else:
                raise

    def _checkversion(self):
        import pytest
        minver = self.inicfg.get('minversion', None)
        if minver:
            ver = minver.split(".")
            myver = pytest.__version__.split(".")
            if myver < ver:
                raise pytest.UsageError(
                    "%s:%d: requires pytest-%s, actual pytest-%s'" % (
                        self.inicfg.config.path, self.inicfg.lineof('minversion'),
                        minver, pytest.__version__))

    def parse(self, args, addopts=True):
        # parse given cmdline arguments into this config object.
        assert not hasattr(self, 'args'), (
            "can only parse cmdline args at most once per Config object")
        self._origargs = args
        self.hook.pytest_addhooks.call_historic(
            kwargs=dict(pluginmanager=self.pluginmanager))
        self._preparse(args, addopts=addopts)
        # XXX deprecated hook:
        self.hook.pytest_cmdline_preparse(config=self, args=args)
        self._parser.after_preparse = True
        try:
            args = self._parser.parse_setoption(args, self.option, namespace=self.option)
            if not args:
                cwd = os.getcwd()
                if cwd == self.rootdir:
                    args = self.getini('testpaths')
                if not args:
                    args = [cwd]
            self.args = args
        except PrintHelp:
            pass

    def addinivalue_line(self, name, line):
        """ add a line to an ini-file option. The option must have been
        declared but might not yet be set in which case the line becomes the
        the first line in its value. """
        x = self.getini(name)
        assert isinstance(x, list)
        x.append(line)  # modifies the cached list inline

    def getini(self, name):
        """ return configuration value from an :ref:`ini file <inifiles>`. If the
        specified name hasn't been registered through a prior
        :py:func:`parser.addini <_pytest.config.Parser.addini>`
        call (usually from a plugin), a ValueError is raised. """
        try:
            return self._inicache[name]
        except KeyError:
            self._inicache[name] = val = self._getini(name)
            return val

    def _getini(self, name):
        try:
            description, type, default = self._parser._inidict[name]
        except KeyError:
            raise ValueError("unknown configuration value: %r" % (name,))
        value = self._get_override_ini_value(name)
        if value is None:
            try:
                value = self.inicfg[name]
            except KeyError:
                if default is not None:
                    return default
                if type is None:
                    return ''
                return []
        if type == "pathlist":
            dp = py.path.local(self.inicfg.config.path).dirpath()
            values = []
            for relpath in shlex.split(value):
                values.append(dp.join(relpath, abs=True))
            return values
        elif type == "args":
            return shlex.split(value)
        elif type == "linelist":
            return [t for t in map(lambda x: x.strip(), value.split("\n")) if t]
        elif type == "bool":
            return bool(_strtobool(value.strip()))
        else:
            assert type is None
            return value

    def _getconftest_pathlist(self, name, path):
        try:
            mod, relroots = self.pluginmanager._rget_with_confmod(name, path)
        except KeyError:
            return None
        modpath = py.path.local(mod.__file__).dirpath()
        values = []
        for relroot in relroots:
            if not isinstance(relroot, py.path.local):
                relroot = relroot.replace("/", py.path.local.sep)
                relroot = modpath.join(relroot, abs=True)
            values.append(relroot)
        return values

    def _get_override_ini_value(self, name):
        value = None
        # override_ini is a list of "ini=value" options
        # always use the last item if multiple values are set for same ini-name,
        # e.g. -o foo=bar1 -o foo=bar2 will set foo to bar2
        for ini_config in self._override_ini:
            try:
                key, user_ini_value = ini_config.split("=", 1)
            except ValueError:
                raise UsageError("-o/--override-ini expects option=value style.")
            else:
                if key == name:
                    value = user_ini_value
        return value

    def getoption(self, name, default=notset, skip=False):
        """ return command line option value.

        :arg name: name of the option.  You may also specify
            the literal ``--OPT`` option instead of the "dest" option name.
        :arg default: default value if no option of that name exists.
        :arg skip: if True raise pytest.skip if option does not exists
            or has a None value.
        """
        name = self._opt2dest.get(name, name)
        try:
            val = getattr(self.option, name)
            if val is None and skip:
                raise AttributeError(name)
            return val
        except AttributeError:
            if default is not notset:
                return default
            if skip:
                import pytest
                pytest.skip("no %r option found" % (name,))
            raise ValueError("no option named %r" % (name,))

    def getvalue(self, name, path=None):
        """ (deprecated, use getoption()) """
        return self.getoption(name)

    def getvalueorskip(self, name, path=None):
        """ (deprecated, use getoption(skip=True)) """
        return self.getoption(name, skip=True)


def _assertion_supported():
    try:
        assert False
    except AssertionError:
        return True
    else:
        return False


def _warn_about_missing_assertion(mode):
    if not _assertion_supported():
        if mode == 'plain':
            sys.stderr.write("WARNING: ASSERTIONS ARE NOT EXECUTED"
                             " and FAILING TESTS WILL PASS.  Are you"
                             " using python -O?")
        else:
            sys.stderr.write("WARNING: assertions not in test modules or"
                             " plugins will be ignored"
                             " because assert statements are not executed "
                             "by the underlying Python interpreter "
                             "(are you using python -O?)\n")


def exists(path, ignore=EnvironmentError):
    try:
        return path.check()
    except ignore:
        return False


def getcfg(args, warnfunc=None):
    """
    Search the list of arguments for a valid ini-file for pytest,
    and return a tuple of (rootdir, inifile, cfg-dict).

    note: warnfunc is an optional function used to warn
        about ini-files that use deprecated features.
        This parameter should be removed when pytest
        adopts standard deprecation warnings (#1804).
    """
    from _pytest.deprecated import CFG_PYTEST_SECTION
    inibasenames = ["pytest.ini", "tox.ini", "setup.cfg"]
    args = [x for x in args if not str(x).startswith("-")]
    if not args:
        args = [py.path.local()]
    for arg in args:
        arg = py.path.local(arg)
        for base in arg.parts(reverse=True):
            for inibasename in inibasenames:
                p = base.join(inibasename)
                if exists(p):
                    iniconfig = py.iniconfig.IniConfig(p)
                    if 'pytest' in iniconfig.sections:
                        if inibasename == 'setup.cfg' and warnfunc:
                            warnfunc('C1', CFG_PYTEST_SECTION.format(filename=inibasename))
                        return base, p, iniconfig['pytest']
                    if inibasename == 'setup.cfg' and 'tool:pytest' in iniconfig.sections:
                        return base, p, iniconfig['tool:pytest']
                    elif inibasename == "pytest.ini":
                        # allowed to be empty
                        return base, p, {}
    return None, None, None


def get_common_ancestor(paths):
    common_ancestor = None
    for path in paths:
        if not path.exists():
            continue
        if common_ancestor is None:
            common_ancestor = path
        else:
            if path.relto(common_ancestor) or path == common_ancestor:
                continue
            elif common_ancestor.relto(path):
                common_ancestor = path
            else:
                shared = path.common(common_ancestor)
                if shared is not None:
                    common_ancestor = shared
    if common_ancestor is None:
        common_ancestor = py.path.local()
    elif common_ancestor.isfile():
        common_ancestor = common_ancestor.dirpath()
    return common_ancestor


def get_dirs_from_args(args):
    def is_option(x):
        return str(x).startswith('-')

    def get_file_part_from_node_id(x):
        return str(x).split('::')[0]

    def get_dir_from_path(path):
        if path.isdir():
            return path
        return py.path.local(path.dirname)

    # These look like paths but may not exist
    possible_paths = (
        py.path.local(get_file_part_from_node_id(arg))
        for arg in args
        if not is_option(arg)
    )

    return [
        get_dir_from_path(path)
        for path in possible_paths
        if path.exists()
    ]


def determine_setup(inifile, args, warnfunc=None, rootdir_cmd_arg=None):
    dirs = get_dirs_from_args(args)
    if inifile:
        iniconfig = py.iniconfig.IniConfig(inifile)
        is_cfg_file = str(inifile).endswith('.cfg')
        # TODO: [pytest] section in *.cfg files is depricated. Need refactoring.
        sections = ['tool:pytest', 'pytest'] if is_cfg_file else ['pytest']
        for section in sections:
            try:
                inicfg = iniconfig[section]
                if is_cfg_file and section == 'pytest' and warnfunc:
                    from _pytest.deprecated import CFG_PYTEST_SECTION
                    warnfunc('C1', CFG_PYTEST_SECTION.format(filename=str(inifile)))
                break
            except KeyError:
                inicfg = None
        rootdir = get_common_ancestor(dirs)
    else:
        ancestor = get_common_ancestor(dirs)
        rootdir, inifile, inicfg = getcfg([ancestor], warnfunc=warnfunc)
        if rootdir is None:
            for rootdir in ancestor.parts(reverse=True):
                if rootdir.join("setup.py").exists():
                    break
            else:
                rootdir, inifile, inicfg = getcfg(dirs, warnfunc=warnfunc)
                if rootdir is None:
                    rootdir = get_common_ancestor([py.path.local(), ancestor])
                    is_fs_root = os.path.splitdrive(str(rootdir))[1] == '/'
                    if is_fs_root:
                        rootdir = ancestor
    if rootdir_cmd_arg:
        rootdir_abs_path = py.path.local(os.path.expandvars(rootdir_cmd_arg))
        if not os.path.isdir(str(rootdir_abs_path)):
            raise UsageError("Directory '{}' not found. Check your '--rootdir' option.".format(rootdir_abs_path))
        rootdir = rootdir_abs_path
    return rootdir, inifile, inicfg or {}


def setns(obj, dic):
    import pytest
    for name, value in dic.items():
        if isinstance(value, dict):
            mod = getattr(obj, name, None)
            if mod is None:
                modname = "pytest.%s" % name
                mod = types.ModuleType(modname)
                sys.modules[modname] = mod
                mod.__all__ = []
                setattr(obj, name, mod)
            obj.__all__.append(name)
            setns(mod, value)
        else:
            setattr(obj, name, value)
            obj.__all__.append(name)
            # if obj != pytest:
            #    pytest.__all__.append(name)
            setattr(pytest, name, value)


def create_terminal_writer(config, *args, **kwargs):
    """Create a TerminalWriter instance configured according to the options
    in the config object. Every code which requires a TerminalWriter object
    and has access to a config object should use this function.
    """
    tw = py.io.TerminalWriter(*args, **kwargs)
    if config.option.color == 'yes':
        tw.hasmarkup = True
    if config.option.color == 'no':
        tw.hasmarkup = False
    return tw


def _strtobool(val):
    """Convert a string representation of truth to true (1) or false (0).

    True values are 'y', 'yes', 't', 'true', 'on', and '1'; false values
    are 'n', 'no', 'f', 'false', 'off', and '0'.  Raises ValueError if
    'val' is anything else.

    .. note:: copied from distutils.util
    """
    val = val.lower()
    if val in ('y', 'yes', 't', 'true', 'on', '1'):
        return 1
    elif val in ('n', 'no', 'f', 'false', 'off', '0'):
        return 0
    else:
        raise ValueError("invalid truth value %r" % (val,))
