from __future__ import absolute_import, division, print_function

import pytest

import matplotlib
from matplotlib.testing.decorators import image_comparison
import matplotlib.pyplot as plt


@pytest.mark.skipif(not matplotlib.checkdep_usetex(True),
                    reason='Missing TeX or Ghostscript or dvipng')
@image_comparison(baseline_images=['test_usetex'],
                  extensions=['pdf', 'png'],
                  tol=0.3)
def test_usetex():
    matplotlib.rcParams['text.usetex'] = True
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.text(0.1, 0.2,
            # the \LaTeX macro exercises character sizing and placement,
            # \left[ ... \right\} draw some variable-height characters,
            # \sqrt and \frac draw horizontal rules, \mathrm changes the font
            r'\LaTeX\ $\left[\int\limits_e^{2e}'
            r'\sqrt\frac{\log^3 x}{x}\,\mathrm{d}x \right\}$',
            fontsize=24)
    ax.set_xticks([])
    ax.set_yticks([])
