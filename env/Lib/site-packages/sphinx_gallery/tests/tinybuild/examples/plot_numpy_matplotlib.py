"""
======================
Link to other packages
======================

Use :mod:`sphinx_gallery` to link to other packages, like
:mod:`numpy`, :mod:`matplotlib.colors`, and :mod:`matplotlib.pyplot`.

FYI this gallery uses :obj:`sphinx_gallery.sorting.FileNameSortKey`.
"""

from itertools import compress  # noqa
from warnings import warn

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from local_module import N  # N = 1000
from matplotlib.colors import is_color_like
from matplotlib.figure import Figure

import sphinx_gallery.backreferences
import sphinx_gallery._dummy.nested
from sphinx_gallery.py_source_parser import Block

t = np.arange(N) / float(N)
win = np.hanning(N)
print(is_color_like("r"))
fig, ax = plt.subplots()
ax.plot(t, win, color="r")
ax.text(0, 1, "png", size=40, va="top")
fig.tight_layout()
orig_dpi = 80.0 if matplotlib.__version__[0] < "2" else 100.0
assert plt.rcParams["figure.dpi"] == orig_dpi
plt.rcParams["figure.dpi"] = 70.0
assert plt.rcParams["figure.dpi"] == 70.0
fig3d, ax3d = plt.subplots(subplot_kw={"projection": "3d"})
ax3d.plot(t, t, win, color="r")  # This should map to Axes3D.plot, not Axes.plot.
listy = [0, 1]
compress("abc", [0, 0, 1])
warn("This warning should show up in the output", RuntimeWarning)
x = Figure()  # plt.Figure should be decorated (class), x shouldn't (inst)
# nested resolution resolves to numpy.random.mtrand.RandomState:
rng = np.random.RandomState(0)
# test Issue 583
sphinx_gallery.backreferences.identify_names(
    [Block("text", "Text block", 1)],
    sphinx_gallery.backreferences._make_ref_regex(),
)
# 583: methods don't link properly
dc = sphinx_gallery._dummy.DummyClass()
dc.run()
print(dc.prop)
# 1364: nested methods don't link properly
ndc = sphinx_gallery._dummy.nested.NestedDummyClass()
ndc.run()
print(ndc.prop)
