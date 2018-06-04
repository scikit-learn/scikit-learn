import importlib
import os
import sys

from matplotlib.compat.subprocess import Popen
import pytest


# Minimal smoke-testing of the backends for which the dependencies are
# PyPI-installable on Travis.  They are not available for all tested Python
# versions so we don't fail on missing backends.
#
# We also don't test on Py2 because its subprocess module doesn't support
# timeouts, and it would require a separate code path to check for module
# existence without actually trying to import the module (which may install
# an undesirable input hook).


def _get_testable_interactive_backends():
    backends = []
    for deps, backend in [(["cairocffi", "pgi"], "gtk3agg"),
                          (["cairocffi", "pgi"], "gtk3cairo"),
                          (["PyQt5"], "qt5agg"),
                          (["cairocffi", "PyQt5"], "qt5cairo"),
                          (["tkinter"], "tkagg"),
                          (["wx"], "wxagg")]:
        reason = None
        if sys.version_info < (3,):
            reason = "Py3-only test"
        elif not os.environ.get("DISPLAY"):
            reason = "No $DISPLAY"
        elif any(importlib.util.find_spec(dep) is None for dep in deps):
            reason = "Missing dependency"
        backends.append(pytest.mark.skip(reason=reason)(backend) if reason
                        else backend)
    return backends


_test_script = """\
import sys
from matplotlib import pyplot as plt

fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot([1,2,3], [1,3,1])
fig.canvas.mpl_connect("draw_event", lambda event: sys.exit())
plt.show()
"""


@pytest.mark.parametrize("backend", _get_testable_interactive_backends())
@pytest.mark.flaky(reruns=3)
def test_backend(backend):
    environ = os.environ.copy()
    environ["MPLBACKEND"] = backend
    proc = Popen([sys.executable, "-c", _test_script], env=environ)
    # Empirically, 1s is not enough on Travis.
    assert proc.wait(timeout=10) == 0
