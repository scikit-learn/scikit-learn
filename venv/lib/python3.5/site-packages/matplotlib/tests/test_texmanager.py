from __future__ import absolute_import, division, print_function

import matplotlib.pyplot as plt
from matplotlib.texmanager import TexManager


def test_fontconfig_preamble():
    """
    Test that the preamble is included in _fontconfig
    """
    plt.rcParams['text.usetex'] = True

    tm1 = TexManager()
    font_config1 = tm1.get_font_config()

    plt.rcParams['text.latex.preamble'] = ['\\usepackage{txfonts}']
    tm2 = TexManager()
    font_config2 = tm2.get_font_config()

    assert font_config1 != font_config2
