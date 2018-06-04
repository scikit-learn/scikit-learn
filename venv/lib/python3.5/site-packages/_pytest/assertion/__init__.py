"""
support for presenting detailed information in failing assertions.
"""
from __future__ import absolute_import, division, print_function
import sys
import six

from _pytest.assertion import util
from _pytest.assertion import rewrite
from _pytest.assertion import truncate


def pytest_addoption(parser):
    group = parser.getgroup("debugconfig")
    group.addoption('--assert',
                    action="store",
                    dest="assertmode",
                    choices=("rewrite", "plain",),
                    default="rewrite",
                    metavar="MODE",
                    help="""Control assertion debugging tools.  'plain'
                            performs no assertion debugging.  'rewrite'
                            (the default) rewrites assert statements in
                            test modules on import to provide assert
                            expression information.""")


def register_assert_rewrite(*names):
    """Register one or more module names to be rewritten on import.

    This function will make sure that this module or all modules inside
    the package will get their assert statements rewritten.
    Thus you should make sure to call this before the module is
    actually imported, usually in your __init__.py if you are a plugin
    using a package.

    :raise TypeError: if the given module names are not strings.
    """
    for name in names:
        if not isinstance(name, str):
            msg = 'expected module names as *args, got {0} instead'
            raise TypeError(msg.format(repr(names)))
    for hook in sys.meta_path:
        if isinstance(hook, rewrite.AssertionRewritingHook):
            importhook = hook
            break
    else:
        importhook = DummyRewriteHook()
    importhook.mark_rewrite(*names)


class DummyRewriteHook(object):
    """A no-op import hook for when rewriting is disabled."""

    def mark_rewrite(self, *names):
        pass


class AssertionState(object):
    """State for the assertion plugin."""

    def __init__(self, config, mode):
        self.mode = mode
        self.trace = config.trace.root.get("assertion")
        self.hook = None


def install_importhook(config):
    """Try to install the rewrite hook, raise SystemError if it fails."""
    # Jython has an AST bug that make the assertion rewriting hook malfunction.
    if (sys.platform.startswith('java')):
        raise SystemError('rewrite not supported')

    config._assertstate = AssertionState(config, 'rewrite')
    config._assertstate.hook = hook = rewrite.AssertionRewritingHook(config)
    sys.meta_path.insert(0, hook)
    config._assertstate.trace('installed rewrite import hook')

    def undo():
        hook = config._assertstate.hook
        if hook is not None and hook in sys.meta_path:
            sys.meta_path.remove(hook)

    config.add_cleanup(undo)
    return hook


def pytest_collection(session):
    # this hook is only called when test modules are collected
    # so for example not in the master process of pytest-xdist
    # (which does not collect test modules)
    assertstate = getattr(session.config, '_assertstate', None)
    if assertstate:
        if assertstate.hook is not None:
            assertstate.hook.set_session(session)


def pytest_runtest_setup(item):
    """Setup the pytest_assertrepr_compare hook

    The newinterpret and rewrite modules will use util._reprcompare if
    it exists to use custom reporting via the
    pytest_assertrepr_compare hook.  This sets up this custom
    comparison for the test.
    """
    def callbinrepr(op, left, right):
        """Call the pytest_assertrepr_compare hook and prepare the result

        This uses the first result from the hook and then ensures the
        following:
        * Overly verbose explanations are truncated unless configured otherwise
          (eg. if running in verbose mode).
        * Embedded newlines are escaped to help util.format_explanation()
          later.
        * If the rewrite mode is used embedded %-characters are replaced
          to protect later % formatting.

        The result can be formatted by util.format_explanation() for
        pretty printing.
        """
        hook_result = item.ihook.pytest_assertrepr_compare(
            config=item.config, op=op, left=left, right=right)
        for new_expl in hook_result:
            if new_expl:
                new_expl = truncate.truncate_if_required(new_expl, item)
                new_expl = [line.replace("\n", "\\n") for line in new_expl]
                res = six.text_type("\n~").join(new_expl)
                if item.config.getvalue("assertmode") == "rewrite":
                    res = res.replace("%", "%%")
                return res
    util._reprcompare = callbinrepr


def pytest_runtest_teardown(item):
    util._reprcompare = None


def pytest_sessionfinish(session):
    assertstate = getattr(session.config, '_assertstate', None)
    if assertstate:
        if assertstate.hook is not None:
            assertstate.hook.set_session(None)


# Expose this plugin's implementation for the pytest_assertrepr_compare hook
pytest_assertrepr_compare = util.assertrepr_compare
