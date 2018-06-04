from __future__ import absolute_import, division, print_function

import six

import os
import tempfile
import warnings

import numpy as np
import pytest

from matplotlib.font_manager import (
    findfont, FontProperties, fontManager, json_dump, json_load, get_font,
    get_fontconfig_fonts, is_opentype_cff_font, fontManager as fm)
from matplotlib import rc_context

if six.PY2:
    from distutils.spawn import find_executable
    has_fclist = find_executable('fc-list') is not None
else:
    # py >= 3.3
    from shutil import which
    has_fclist = which('fc-list') is not None


def test_font_priority():
    with rc_context(rc={
            'font.sans-serif':
            ['cmmi10', 'Bitstream Vera Sans']}):
        font = findfont(
            FontProperties(family=["sans-serif"]))
    assert os.path.basename(font) == 'cmmi10.ttf'

    # Smoketest get_charmap, which isn't used internally anymore
    font = get_font(font)
    cmap = font.get_charmap()
    assert len(cmap) == 131
    assert cmap[8729] == 30


def test_score_weight():
    assert 0 == fontManager.score_weight("regular", "regular")
    assert 0 == fontManager.score_weight("bold", "bold")
    assert (0 < fontManager.score_weight(400, 400) <
            fontManager.score_weight("normal", "bold"))
    assert (0 < fontManager.score_weight("normal", "regular") <
            fontManager.score_weight("normal", "bold"))
    assert (fontManager.score_weight("normal", "regular") ==
            fontManager.score_weight(400, 400))


def test_json_serialization():
    # on windows, we can't open a file twice, so save the name and unlink
    # manually...
    try:
        name = None
        with tempfile.NamedTemporaryFile(delete=False) as temp:
            name = temp.name
        json_dump(fontManager, name)
        copy = json_load(name)
    finally:
        if name and os.path.exists(name):
            os.remove(name)
    with warnings.catch_warnings():
        warnings.filterwarnings('ignore', 'findfont: Font family.*not found')
        for prop in ({'family': 'STIXGeneral'},
                     {'family': 'Bitstream Vera Sans', 'weight': 700},
                     {'family': 'no such font family'}):
            fp = FontProperties(**prop)
            assert (fontManager.findfont(fp, rebuild_if_missing=False) ==
                    copy.findfont(fp, rebuild_if_missing=False))


def test_otf():
    fname = '/usr/share/fonts/opentype/freefont/FreeMono.otf'
    if os.path.exists(fname):
        assert is_opentype_cff_font(fname)

    otf_files = [f for f in fm.ttffiles if 'otf' in f]
    for f in otf_files:
        with open(f, 'rb') as fd:
            res = fd.read(4) == b'OTTO'
        assert res == is_opentype_cff_font(f)


@pytest.mark.skipif(not has_fclist, reason='no fontconfig installed')
def test_get_fontconfig_fonts():
    assert len(get_fontconfig_fonts()) > 1


@pytest.mark.parametrize('factor', [2, 4, 6, 8])
def test_hinting_factor(factor):
    font = findfont(FontProperties(family=["sans-serif"]))

    font1 = get_font(font, hinting_factor=1)
    font1.clear()
    font1.set_size(12, 100)
    font1.set_text('abc')
    expected = font1.get_width_height()

    hinted_font = get_font(font, hinting_factor=factor)
    hinted_font.clear()
    hinted_font.set_size(12, 100)
    hinted_font.set_text('abc')
    # Check that hinting only changes text layout by a small (10%) amount.
    np.testing.assert_allclose(hinted_font.get_width_height(), expected,
                               rtol=0.1)
