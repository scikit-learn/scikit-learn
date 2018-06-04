# -*- coding: utf-8 -*-

from __future__ import absolute_import, division, print_function
from six import BytesIO

import matplotlib.afm as afm


AFM_TEST_DATA = b"""StartFontMetrics 2.0
Comment Comments are ignored.
Comment Creation Date:Mon Nov 13 12:34:11 GMT 2017
FontName MyFont-Bold
EncodingScheme FontSpecific
FullName My Font Bold
FamilyName Test Fonts
Weight Bold
ItalicAngle 0.0
IsFixedPitch false
UnderlinePosition -100
UnderlineThickness 50
Version 001.000
Notice Copyright (c) 2017 No one.
FontBBox 0 -321 1234 369
StartCharMetrics 3
C 0 ; WX 250 ; N space ; B 0 0 0 0 ;
C 42 ; WX 1141 ; N foo ; B 40 60 800 360 ;
C 99 ; WX 583 ; N bar ; B 40 -10 543 210 ;
EndCharMetrics
EndFontMetrics
"""


def test_nonascii_str():
    # This tests that we also decode bytes as utf-8 properly.
    # Else, font files with non ascii characters fail to load.
    inp_str = u"привет"
    byte_str = inp_str.encode("utf8")

    ret = afm._to_str(byte_str)
    assert ret == inp_str


def test_parse_header():
    fh = BytesIO(AFM_TEST_DATA)
    header = afm._parse_header(fh)
    assert header == {
        b'StartFontMetrics': 2.0,
        b'FontName': 'MyFont-Bold',
        b'EncodingScheme': 'FontSpecific',
        b'FullName': 'My Font Bold',
        b'FamilyName': 'Test Fonts',
        b'Weight': 'Bold',
        b'ItalicAngle': 0.0,
        b'IsFixedPitch': False,
        b'UnderlinePosition': -100,
        b'UnderlineThickness': 50,
        b'Version': '001.000',
        b'Notice': 'Copyright (c) 2017 No one.',
        b'FontBBox': [0, -321, 1234, 369],
        b'StartCharMetrics': 3,
    }


def test_parse_char_metrics():
    fh = BytesIO(AFM_TEST_DATA)
    afm._parse_header(fh)  # position
    metrics = afm._parse_char_metrics(fh)
    assert metrics == (
        {0: (250.0, 'space', [0, 0, 0, 0]),
         42: (1141.0, 'foo', [40, 60, 800, 360]),
         99: (583.0, 'bar', [40, -10, 543, 210]),
         },
        {'space': (250.0, [0, 0, 0, 0]),
         'foo': (1141.0, [40, 60, 800, 360]),
         'bar': (583.0, [40, -10, 543, 210]),
         })


def test_get_familyname_guessed():
    fh = BytesIO(AFM_TEST_DATA)
    fm = afm.AFM(fh)
    del fm._header[b'FamilyName']  # remove FamilyName, so we have to guess
    assert fm.get_familyname() == 'My Font'
