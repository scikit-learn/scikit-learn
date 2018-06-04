from __future__ import absolute_import, division, print_function

import six
from matplotlib.testing.decorators import skip_if_command_unavailable

import matplotlib.dviread as dr
import os.path
import json
import pytest


def test_PsfontsMap(monkeypatch):
    monkeypatch.setattr(dr, 'find_tex_file', lambda x: x)

    filename = os.path.join(
        os.path.dirname(__file__),
        'baseline_images', 'dviread', 'test.map')
    fontmap = dr.PsfontsMap(filename)
    # Check all properties of a few fonts
    for n in [1, 2, 3, 4, 5]:
        key = ('TeXfont%d' % n).encode('ascii')
        entry = fontmap[key]
        assert entry.texname == key
        assert entry.psname == ('PSfont%d' % n).encode('ascii')
        if n not in [3, 5]:
            assert entry.encoding == ('font%d.enc' % n).encode('ascii')
        elif n == 3:
            assert entry.encoding == b'enc3.foo'
        # We don't care about the encoding of TeXfont5, which specifies
        # multiple encodings.
        if n not in [1, 5]:
            assert entry.filename == ('font%d.pfa' % n).encode('ascii')
        else:
            assert entry.filename == ('font%d.pfb' % n).encode('ascii')
        if n == 4:
            assert entry.effects == {'slant': -0.1, 'extend': 2.2}
        else:
            assert entry.effects == {}
    # Some special cases
    entry = fontmap[b'TeXfont6']
    assert entry.filename is None
    assert entry.encoding is None
    entry = fontmap[b'TeXfont7']
    assert entry.filename is None
    assert entry.encoding == b'font7.enc'
    entry = fontmap[b'TeXfont8']
    assert entry.filename == b'font8.pfb'
    assert entry.encoding is None
    entry = fontmap[b'TeXfont9']
    assert entry.filename == b'/absolute/font9.pfb'
    # Missing font
    with pytest.raises(KeyError) as exc:
        fontmap[b'no-such-font']
    assert 'no-such-font' in str(exc.value)


@skip_if_command_unavailable(["kpsewhich", "-version"])
def test_dviread():
    dir = os.path.join(os.path.dirname(__file__), 'baseline_images', 'dviread')
    with open(os.path.join(dir, 'test.json')) as f:
        correct = json.load(f)
        for entry in correct:
            entry['text'] = [[a, b, c, d.encode('ascii'), e]
                             for [a, b, c, d, e] in entry['text']]
    with dr.Dvi(os.path.join(dir, 'test.dvi'), None) as dvi:
        data = [{'text': [[t.x, t.y,
                           six.unichr(t.glyph),
                           t.font.texname,
                           round(t.font.size, 2)]
                          for t in page.text],
                 'boxes': [[b.x, b.y, b.height, b.width] for b in page.boxes]}
                for page in dvi]
    assert data == correct
