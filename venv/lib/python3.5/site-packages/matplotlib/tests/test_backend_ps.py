# -*- coding: utf-8 -*-

from __future__ import absolute_import, division, print_function

import io
import re

import numpy as np
import pytest
import six

import matplotlib
import matplotlib.pyplot as plt
from matplotlib import patheffects
from matplotlib.testing.decorators import image_comparison
from matplotlib.testing.determinism import (_determinism_source_date_epoch,
                                            _determinism_check)


needs_ghostscript = pytest.mark.xfail(
    matplotlib.checkdep_ghostscript()[0] is None,
    reason="This test needs a ghostscript installation")


needs_usetex = pytest.mark.xfail(
    not matplotlib.checkdep_usetex(True),
    reason="This test needs a TeX installation")


# This tests tends to hit a TeX cache lock on AppVeyor.
@pytest.mark.flaky(reruns=3)
@pytest.mark.parametrize('format, use_log, rcParams', [
    ('ps', False, {}),
    needs_ghostscript(('ps', False, {'ps.usedistiller': 'ghostscript'})),
    needs_usetex(needs_ghostscript(('ps', False, {'text.latex.unicode': True,
                                                  'text.usetex': True}))),
    ('eps', False, {}),
    ('eps', True, {'ps.useafm': True}),
    needs_usetex(needs_ghostscript(('eps', False, {'text.latex.unicode': True,
                                                   'text.usetex': True}))),
], ids=[
    'ps',
    'ps with distiller',
    'ps with usetex',
    'eps',
    'eps afm',
    'eps with usetex'
])
def test_savefig_to_stringio(format, use_log, rcParams):
    matplotlib.rcParams.update(rcParams)

    fig, ax = plt.subplots()
    buffers = [
        six.moves.StringIO(),
        io.StringIO(),
        io.BytesIO()]

    if use_log:
        ax.set_yscale('log')

    ax.plot([1, 2], [1, 2])
    ax.set_title(u"Déjà vu")
    for buffer in buffers:
        fig.savefig(buffer, format=format)

    values = [x.getvalue() for x in buffers]

    if six.PY3:
        values = [
            values[0].encode('ascii'),
            values[1].encode('ascii'),
            values[2]]

    # Remove comments from the output.  This includes things that
    # could change from run to run, such as the time.
    values = [re.sub(b'%%.*?\n', b'', x) for x in values]

    assert values[0] == values[1]
    assert values[1] == values[2].replace(b'\r\n', b'\n')
    for buffer in buffers:
        buffer.close()


def test_patheffects():
    with matplotlib.rc_context():
        matplotlib.rcParams['path.effects'] = [
            patheffects.withStroke(linewidth=4, foreground='w')]
        fig, ax = plt.subplots()
        ax.plot([1, 2, 3])
        with io.BytesIO() as ps:
            fig.savefig(ps, format='ps')


@needs_usetex
@needs_ghostscript
def test_tilde_in_tempfilename():
    # Tilde ~ in the tempdir path (e.g. TMPDIR, TMP or TEMP on windows
    # when the username is very long and windows uses a short name) breaks
    # latex before https://github.com/matplotlib/matplotlib/pull/5928
    import tempfile
    import shutil
    import os
    import os.path

    tempdir = None
    old_tempdir = tempfile.tempdir
    try:
        # change the path for new tempdirs, which is used
        # internally by the ps backend to write a file
        tempdir = tempfile.mkdtemp()
        base_tempdir = os.path.join(tempdir, "short~1")
        os.makedirs(base_tempdir)
        tempfile.tempdir = base_tempdir

        # usetex results in the latex call, which does not like the ~
        plt.rc('text', usetex=True)
        plt.plot([1, 2, 3, 4])
        plt.xlabel(r'\textbf{time} (s)')
        output_eps = os.path.join(base_tempdir, 'tex_demo.eps')
        # use the PS backend to write the file...
        plt.savefig(output_eps, format="ps")
    finally:
        tempfile.tempdir = old_tempdir
        if tempdir:
            try:
                shutil.rmtree(tempdir)
            except Exception as e:
                # do not break if this is not removable...
                print(e)


def test_source_date_epoch():
    """Test SOURCE_DATE_EPOCH support for PS output"""
    # SOURCE_DATE_EPOCH support is not tested with text.usetex,
    # because the produced timestamp comes from ghostscript:
    # %%CreationDate: D:20000101000000Z00\'00\', and this could change
    # with another ghostscript version.
    _determinism_source_date_epoch(
        "ps", b"%%CreationDate: Sat Jan 01 00:00:00 2000")


def test_determinism_all():
    """Test for reproducible PS output"""
    _determinism_check(format="ps")


@needs_usetex
@needs_ghostscript
def test_determinism_all_tex():
    """Test for reproducible PS/tex output"""
    _determinism_check(format="ps", usetex=True)


@image_comparison(baseline_images=["empty"], extensions=["eps"])
def test_transparency():
    fig, ax = plt.subplots()
    ax.set_axis_off()
    ax.plot([0, 1], color="r", alpha=0)
    ax.text(.5, .5, "foo", color="r", alpha=0)


@needs_usetex
def test_failing_latex(tmpdir):
    """Test failing latex subprocess call"""
    path = str(tmpdir.join("tmpoutput.ps"))

    matplotlib.rcParams['text.usetex'] = True

    # This fails with "Double subscript"
    plt.xlabel("$22_2_2$")
    with pytest.raises(RuntimeError):
        plt.savefig(path)
