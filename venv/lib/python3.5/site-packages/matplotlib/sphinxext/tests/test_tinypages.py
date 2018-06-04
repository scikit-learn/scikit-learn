""" Tests for tinypages build using sphinx extensions """

import filecmp
from os.path import join as pjoin, dirname, isdir
from subprocess import call, Popen, PIPE
import sys

import pytest

from matplotlib import cbook


needs_sphinx = pytest.mark.skipif(
    call([sys.executable, '-msphinx', '--help'], stdout=PIPE, stderr=PIPE),
    reason="'{} -msphinx' does not return 0".format(sys.executable))


@cbook.deprecated("2.1", alternative="filecmp.cmp")
def file_same(file1, file2):
    with open(file1, 'rb') as fobj:
        contents1 = fobj.read()
    with open(file2, 'rb') as fobj:
        contents2 = fobj.read()
    return contents1 == contents2


def test_tinypages(tmpdir):
    html_dir = pjoin(str(tmpdir), 'html')
    doctree_dir = pjoin(str(tmpdir), 'doctrees')
    # Build the pages with warnings turned into errors
    cmd = [sys.executable, '-msphinx', '-W', '-b', 'html', '-d', doctree_dir,
           pjoin(dirname(__file__), 'tinypages'), html_dir]
    proc = Popen(cmd, stdout=PIPE, stderr=PIPE)
    out, err = proc.communicate()
    assert proc.returncode == 0, \
        "'{} -msphinx' failed with stdout:\n{}\nstderr:\n{}\n".format(
            sys.executable, out, err)

    assert isdir(html_dir)

    def plot_file(num):
        return pjoin(html_dir, 'some_plots-{0}.png'.format(num))

    range_10, range_6, range_4 = [plot_file(i) for i in range(1, 4)]
    # Plot 5 is range(6) plot
    assert filecmp.cmp(range_6, plot_file(5))
    # Plot 7 is range(4) plot
    assert filecmp.cmp(range_4, plot_file(7))
    # Plot 11 is range(10) plot
    assert filecmp.cmp(range_10, plot_file(11))
    # Plot 12 uses the old range(10) figure and the new range(6) figure
    assert filecmp.cmp(range_10, plot_file('12_00'))
    assert filecmp.cmp(range_6, plot_file('12_01'))
    # Plot 13 shows close-figs in action
    assert filecmp.cmp(range_4, plot_file(13))
    # Plot 14 has included source
    with open(pjoin(html_dir, 'some_plots.html'), 'rb') as fobj:
        html_contents = fobj.read()
    assert b'# Only a comment' in html_contents
    # check plot defined in external file.
    assert filecmp.cmp(range_4, pjoin(html_dir, 'range4.png'))
    assert filecmp.cmp(range_6, pjoin(html_dir, 'range6.png'))
    # check if figure caption made it into html file
    assert b'This is the caption for plot 15.' in html_contents
