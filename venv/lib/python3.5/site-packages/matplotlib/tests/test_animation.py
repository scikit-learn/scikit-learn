from __future__ import absolute_import, division, print_function

import six

import sys
import tempfile

import numpy as np
import pytest

import matplotlib as mpl
from matplotlib import pyplot as plt
from matplotlib import animation


class NullMovieWriter(animation.AbstractMovieWriter):
    """
    A minimal MovieWriter.  It doesn't actually write anything.
    It just saves the arguments that were given to the setup() and
    grab_frame() methods as attributes, and counts how many times
    grab_frame() is called.

    This class doesn't have an __init__ method with the appropriate
    signature, and it doesn't define an isAvailable() method, so
    it cannot be added to the 'writers' registry.
    """

    frame_size_can_vary = True

    def setup(self, fig, outfile, dpi, *args):
        self.fig = fig
        self.outfile = outfile
        self.dpi = dpi
        self.args = args
        self._count = 0

    def grab_frame(self, **savefig_kwargs):
        self.savefig_kwargs = savefig_kwargs
        self._count += 1

    def finish(self):
        pass


def test_null_movie_writer():
    # Test running an animation with NullMovieWriter.

    fig = plt.figure()

    def init():
        pass

    def animate(i):
        pass

    num_frames = 5
    filename = "unused.null"
    dpi = 50
    savefig_kwargs = dict(foo=0)

    anim = animation.FuncAnimation(fig, animate, init_func=init,
                                   frames=num_frames)
    writer = NullMovieWriter()
    anim.save(filename, dpi=dpi, writer=writer,
              savefig_kwargs=savefig_kwargs)

    assert writer.fig == fig
    assert writer.outfile == filename
    assert writer.dpi == dpi
    assert writer.args == ()
    assert writer.savefig_kwargs == savefig_kwargs
    assert writer._count == num_frames


def test_movie_writer_dpi_default():
    # Test setting up movie writer with figure.dpi default.

    fig = plt.figure()

    filename = "unused.null"
    fps = 5
    codec = "unused"
    bitrate = 1
    extra_args = ["unused"]

    def run():
        pass

    writer = animation.MovieWriter(fps, codec, bitrate, extra_args)
    writer._run = run
    writer.setup(fig, filename)
    assert writer.dpi == fig.dpi


@animation.writers.register('null')
class RegisteredNullMovieWriter(NullMovieWriter):

    # To be able to add NullMovieWriter to the 'writers' registry,
    # we must define an __init__ method with a specific signature,
    # and we must define the class method isAvailable().
    # (These methods are not actually required to use an instance
    # of this class as the 'writer' argument of Animation.save().)

    def __init__(self, fps=None, codec=None, bitrate=None,
                 extra_args=None, metadata=None):
        pass

    @classmethod
    def isAvailable(self):
        return True


WRITER_OUTPUT = [
    ('ffmpeg', 'movie.mp4'),
    ('ffmpeg_file', 'movie.mp4'),
    ('avconv', 'movie.mp4'),
    ('avconv_file', 'movie.mp4'),
    ('imagemagick', 'movie.gif'),
    ('imagemagick_file', 'movie.gif'),
    ('pillow', 'movie.gif'),
    ('html', 'movie.html'),
    ('null', 'movie.null')
]
if sys.version_info >= (3, 6):
    from pathlib import Path
    WRITER_OUTPUT += [
        (writer, Path(output)) for writer, output in WRITER_OUTPUT]


# Smoke test for saving animations.  In the future, we should probably
# design more sophisticated tests which compare resulting frames a-la
# matplotlib.testing.image_comparison
@pytest.mark.parametrize('writer, output', WRITER_OUTPUT)
def test_save_animation_smoketest(tmpdir, writer, output):
    try:
        # for ImageMagick the rcparams must be patched to account for
        # 'convert' being a built in MS tool, not the imagemagick
        # tool.
        writer._init_from_registry()
    except AttributeError:
        pass
    if not animation.writers.is_available(writer):
        pytest.skip("writer '%s' not available on this system" % writer)
    fig, ax = plt.subplots()
    line, = ax.plot([], [])

    ax.set_xlim(0, 10)
    ax.set_ylim(-1, 1)

    dpi = None
    codec = None
    if writer == 'ffmpeg':
        # Issue #8253
        fig.set_size_inches((10.85, 9.21))
        dpi = 100.
        codec = 'h264'

    def init():
        line.set_data([], [])
        return line,

    def animate(i):
        x = np.linspace(0, 10, 100)
        y = np.sin(x + i)
        line.set_data(x, y)
        return line,

    # Use temporary directory for the file-based writers, which produce a file
    # per frame with known names.
    with tmpdir.as_cwd():
        anim = animation.FuncAnimation(fig, animate, init_func=init, frames=5)
        try:
            anim.save(output, fps=30, writer=writer, bitrate=500, dpi=dpi,
                      codec=codec)
        except UnicodeDecodeError:
            pytest.xfail("There can be errors in the numpy import stack, "
                         "see issues #1891 and #2679")


def test_no_length_frames():
    fig, ax = plt.subplots()
    line, = ax.plot([], [])

    def init():
        line.set_data([], [])
        return line,

    def animate(i):
        x = np.linspace(0, 10, 100)
        y = np.sin(x + i)
        line.set_data(x, y)
        return line,

    anim = animation.FuncAnimation(fig, animate, init_func=init,
                                   frames=iter(range(5)))
    writer = NullMovieWriter()
    anim.save('unused.null', writer=writer)


def test_movie_writer_registry():
    ffmpeg_path = mpl.rcParams['animation.ffmpeg_path']
    # Not sure about the first state as there could be some writer
    # which set rcparams
    # assert not animation.writers._dirty
    assert len(animation.writers._registered) > 0
    animation.writers.list()  # resets dirty state
    assert not animation.writers._dirty
    mpl.rcParams['animation.ffmpeg_path'] = u"not_available_ever_xxxx"
    assert animation.writers._dirty
    animation.writers.list()  # resets
    assert not animation.writers._dirty
    assert not animation.writers.is_available("ffmpeg")
    # something which is guaranteed to be available in path
    # and exits immediately
    bin = u"true" if sys.platform != 'win32' else u"where"
    mpl.rcParams['animation.ffmpeg_path'] = bin
    assert animation.writers._dirty
    animation.writers.list()  # resets
    assert not animation.writers._dirty
    assert animation.writers.is_available("ffmpeg")
    mpl.rcParams['animation.ffmpeg_path'] = ffmpeg_path
