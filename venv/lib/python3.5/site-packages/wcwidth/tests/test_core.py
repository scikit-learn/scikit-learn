# coding: utf-8
"""Core tests module for wcwidth."""
import wcwidth


def test_hello_jp():
    u"""
    Width of Japanese phrase: コンニチハ, セカイ!

    Given a phrase of 5 and 3 Katakana ideographs, joined with
    3 English-ASCII punctuation characters, totaling 11, this
    phrase consumes 19 cells of a terminal emulator.
    """
    # given,
    phrase = u'コンニチハ, セカイ!'
    expect_length_each = (2, 2, 2, 2, 2, 1, 1, 2, 2, 2, 1)
    expect_length_phrase = sum(expect_length_each)

    # exercise,
    length_each = tuple(map(wcwidth.wcwidth, phrase))
    length_phrase = wcwidth.wcswidth(phrase)

    # verify,
    assert length_each == expect_length_each
    assert length_phrase == expect_length_phrase


def test_wcswidth_substr():
    """
    Test wcswidth() optional 2nd parameter, ``n``.

    ``n`` determines at which position of the string
    to stop counting length.
    """
    # given,
    phrase = u'コンニチハ, セカイ!'
    end = 7
    expect_length_each = (2, 2, 2, 2, 2, 1, 1,)
    expect_length_phrase = sum(expect_length_each)

    # exercise,
    length_phrase = wcwidth.wcswidth(phrase, end)

    # verify,
    assert length_phrase == expect_length_phrase


def test_null_width_0():
    """NULL (0) reports width 0."""
    # given,
    phrase = u'abc\x00def'
    expect_length_each = (1, 1, 1, 0, 1, 1, 1)
    expect_length_phrase = sum(expect_length_each)

    # exercise,
    length_each = tuple(map(wcwidth.wcwidth, phrase))
    length_phrase = wcwidth.wcswidth(phrase, len(phrase))

    # verify,
    assert length_each == expect_length_each
    assert length_phrase == expect_length_phrase


def test_control_c0_width_negative_1():
    """CSI (Control sequence initiate) reports width -1."""
    # given,
    phrase = u'\x1b[0m'
    expect_length_each = (-1, 1, 1, 1)
    expect_length_phrase = -1

    # exercise,
    length_each = tuple(map(wcwidth.wcwidth, phrase))
    length_phrase = wcwidth.wcswidth(phrase, len(phrase))

    # verify,
    assert length_each == expect_length_each
    assert length_phrase == expect_length_phrase


def test_combining_width_negative_1():
    """Simple test combining reports total width of 4."""
    # given,
    phrase = u'--\u05bf--'
    expect_length_each = (1, 1, 0, 1, 1)
    expect_length_phrase = 4

    # exercise,
    length_each = tuple(map(wcwidth.wcwidth, phrase))
    length_phrase = wcwidth.wcswidth(phrase, len(phrase))

    # verify,
    assert length_each == expect_length_each
    assert length_phrase == expect_length_phrase


def test_combining_cafe():
    u"""Phrase cafe + COMBINING ACUTE ACCENT is café of length 4."""
    phrase = u"cafe\u0301"
    expect_length_each = (1, 1, 1, 1, 0)
    expect_length_phrase = 4

    # exercise,
    length_each = tuple(map(wcwidth.wcwidth, phrase))
    length_phrase = wcwidth.wcswidth(phrase, len(phrase))

    # verify,
    assert length_each == expect_length_each
    assert length_phrase == expect_length_phrase


def test_combining_enclosing():
    u"""CYRILLIC CAPITAL LETTER A + COMBINING CYRILLIC HUNDRED THOUSANDS SIGN is А҈ of length 1."""
    phrase = u"\u0410\u0488"
    expect_length_each = (1, 0)
    expect_length_phrase = 1

    # exercise,
    length_each = tuple(map(wcwidth.wcwidth, phrase))
    length_phrase = wcwidth.wcswidth(phrase, len(phrase))

    # verify,
    assert length_each == expect_length_each
    assert length_phrase == expect_length_phrase


def test_combining_spacing():
    u"""Balinese kapal (ship) is ᬓᬨᬮ᭄ of length 4."""
    phrase = u"\u1B13\u1B28\u1B2E\u1B44"
    expect_length_each = (1, 1, 1, 1)
    expect_length_phrase = 4

    # exercise,
    length_each = tuple(map(wcwidth.wcwidth, phrase))
    length_phrase = wcwidth.wcswidth(phrase, len(phrase))

    # verify,
    assert length_each == expect_length_each
    assert length_phrase == expect_length_phrase
