""" version info, help messages, tracing configuration.  """
from __future__ import absolute_import, division, print_function

import py
import pytest
from _pytest.config import PrintHelp
import os
import sys
from argparse import Action


class HelpAction(Action):
    """This is an argparse Action that will raise an exception in
    order to skip the rest of the argument parsing when --help is passed.
    This prevents argparse from quitting due to missing required arguments
    when any are defined, for example by ``pytest_addoption``.
    This is similar to the way that the builtin argparse --help option is
    implemented by raising SystemExit.
    """

    def __init__(self,
                 option_strings,
                 dest=None,
                 default=False,
                 help=None):
        super(HelpAction, self).__init__(
            option_strings=option_strings,
            dest=dest,
            const=True,
            default=default,
            nargs=0,
            help=help)

    def __call__(self, parser, namespace, values, option_string=None):
        setattr(namespace, self.dest, self.const)

        # We should only skip the rest of the parsing after preparse is done
        if getattr(parser._parser, 'after_preparse', False):
            raise PrintHelp


def pytest_addoption(parser):
    group = parser.getgroup('debugconfig')
    group.addoption('--version', action="store_true",
                    help="display pytest lib version and import information.")
    group._addoption("-h", "--help", action=HelpAction, dest="help",
                     help="show help message and configuration info")
    group._addoption('-p', action="append", dest="plugins", default=[],
                     metavar="name",
                     help="early-load given plugin (multi-allowed). "
                     "To avoid loading of plugins, use the `no:` prefix, e.g. "
                     "`no:doctest`.")
    group.addoption('--traceconfig', '--trace-config',
                    action="store_true", default=False,
                    help="trace considerations of conftest.py files."),
    group.addoption('--debug',
                    action="store_true", dest="debug", default=False,
                    help="store internal tracing debug information in 'pytestdebug.log'.")
    group._addoption(
        '-o', '--override-ini', dest="override_ini",
        action="append",
        help='override ini option with "option=value" style, e.g. `-o xfail_strict=True -o cache_dir=cache`.')


@pytest.hookimpl(hookwrapper=True)
def pytest_cmdline_parse():
    outcome = yield
    config = outcome.get_result()
    if config.option.debug:
        path = os.path.abspath("pytestdebug.log")
        debugfile = open(path, 'w')
        debugfile.write("versions pytest-%s, py-%s, "
                        "python-%s\ncwd=%s\nargs=%s\n\n" % (
                            pytest.__version__, py.__version__,
                            ".".join(map(str, sys.version_info)),
                            os.getcwd(), config._origargs))
        config.trace.root.setwriter(debugfile.write)
        undo_tracing = config.pluginmanager.enable_tracing()
        sys.stderr.write("writing pytestdebug information to %s\n" % path)

        def unset_tracing():
            debugfile.close()
            sys.stderr.write("wrote pytestdebug information to %s\n" %
                             debugfile.name)
            config.trace.root.setwriter(None)
            undo_tracing()

        config.add_cleanup(unset_tracing)


def pytest_cmdline_main(config):
    if config.option.version:
        p = py.path.local(pytest.__file__)
        sys.stderr.write("This is pytest version %s, imported from %s\n" %
                         (pytest.__version__, p))
        plugininfo = getpluginversioninfo(config)
        if plugininfo:
            for line in plugininfo:
                sys.stderr.write(line + "\n")
        return 0
    elif config.option.help:
        config._do_configure()
        showhelp(config)
        config._ensure_unconfigure()
        return 0


def showhelp(config):
    reporter = config.pluginmanager.get_plugin('terminalreporter')
    tw = reporter._tw
    tw.write(config._parser.optparser.format_help())
    tw.line()
    tw.line()
    tw.line("[pytest] ini-options in the first "
            "pytest.ini|tox.ini|setup.cfg file found:")
    tw.line()

    for name in config._parser._ininames:
        help, type, default = config._parser._inidict[name]
        if type is None:
            type = "string"
        spec = "%s (%s)" % (name, type)
        line = "  %-24s %s" % (spec, help)
        tw.line(line[:tw.fullwidth])

    tw.line()
    tw.line("environment variables:")
    vars = [
        ("PYTEST_ADDOPTS", "extra command line options"),
        ("PYTEST_PLUGINS", "comma-separated plugins to load during startup"),
        ("PYTEST_DEBUG", "set to enable debug tracing of pytest's internals")
    ]
    for name, help in vars:
        tw.line("  %-24s %s" % (name, help))
    tw.line()
    tw.line()

    tw.line("to see available markers type: pytest --markers")
    tw.line("to see available fixtures type: pytest --fixtures")
    tw.line("(shown according to specified file_or_dir or current dir "
            "if not specified; fixtures with leading '_' are only shown "
            "with the '-v' option")

    for warningreport in reporter.stats.get('warnings', []):
        tw.line("warning : " + warningreport.message, red=True)
    return


conftest_options = [
    ('pytest_plugins', 'list of plugin names to load'),
]


def getpluginversioninfo(config):
    lines = []
    plugininfo = config.pluginmanager.list_plugin_distinfo()
    if plugininfo:
        lines.append("setuptools registered plugins:")
        for plugin, dist in plugininfo:
            loc = getattr(plugin, '__file__', repr(plugin))
            content = "%s-%s at %s" % (dist.project_name, dist.version, loc)
            lines.append("  " + content)
    return lines


def pytest_report_header(config):
    lines = []
    if config.option.debug or config.option.traceconfig:
        lines.append("using: pytest-%s pylib-%s" %
                     (pytest.__version__, py.__version__))

        verinfo = getpluginversioninfo(config)
        if verinfo:
            lines.extend(verinfo)

    if config.option.traceconfig:
        lines.append("active plugins:")
        items = config.pluginmanager.list_name_plugin()
        for name, plugin in items:
            if hasattr(plugin, '__file__'):
                r = plugin.__file__
            else:
                r = repr(plugin)
            lines.append("    %-20s: %s" % (name, r))
    return lines
