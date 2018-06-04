# TODO:
# * Loop Delay is broken on GTKAgg. This is because source_remove() is not
#   working as we want. PyGTK bug?
# * Documentation -- this will need a new section of the User's Guide.
#      Both for Animations and just timers.
#   - Also need to update http://www.scipy.org/Cookbook/Matplotlib/Animations
# * Blit
#   * Currently broken with Qt4 for widgets that don't start on screen
#   * Still a few edge cases that aren't working correctly
#   * Can this integrate better with existing matplotlib animation artist flag?
#     - If animated removes from default draw(), perhaps we could use this to
#       simplify initial draw.
# * Example
#   * Frameless animation - pure procedural with no loop
#   * Need example that uses something like inotify or subprocess
#   * Complex syncing examples
# * Movies
#   * Can blit be enabled for movies?
# * Need to consider event sources to allow clicking through multiple figures
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import six
from six.moves import xrange, zip

import abc
import contextlib
from io import BytesIO
import itertools
import logging
import os
import platform
import sys
import tempfile
import uuid

import numpy as np

from matplotlib._animation_data import (DISPLAY_TEMPLATE, INCLUDED_FRAMES,
                                        JS_INCLUDE)
from matplotlib.compat import subprocess
from matplotlib import cbook, rcParams, rcParamsDefault, rc_context

if six.PY2:
    from base64 import encodestring as encodebytes
else:
    from base64 import encodebytes


_log = logging.getLogger(__name__)

# Process creation flag for subprocess to prevent it raising a terminal
# window. See for example:
# https://stackoverflow.com/questions/24130623/using-python-subprocess-popen-cant-prevent-exe-stopped-working-prompt
if platform.system() == 'Windows':
    subprocess_creation_flags = CREATE_NO_WINDOW = 0x08000000
else:
    # Apparently None won't work here
    subprocess_creation_flags = 0

# Other potential writing methods:
# * http://pymedia.org/
# * libmng (produces swf) python wrappers: https://github.com/libming/libming
# * Wrap x264 API:

# (http://stackoverflow.com/questions/2940671/
# how-to-encode-series-of-images-into-h264-using-x264-api-c-c )


def adjusted_figsize(w, h, dpi, n):
    '''Compute figure size so that pixels are a multiple of n

    Parameters
    ----------
    w, h : float
        Size in inches

    dpi : float
        The dpi

    n : int
        The target multiple

    Returns
    -------
    wnew, hnew : float
        The new figure size in inches.
    '''

    # this maybe simplified if / when we adopt consistent rounding for
    # pixel size across the whole library
    def correct_roundoff(x, dpi, n):
        if int(x*dpi) % n != 0:
            if int(np.nextafter(x, np.inf)*dpi) % n == 0:
                x = np.nextafter(x, np.inf)
            elif int(np.nextafter(x, -np.inf)*dpi) % n == 0:
                x = np.nextafter(x, -np.inf)
        return x

    wnew = int(w * dpi / n) * n / dpi
    hnew = int(h * dpi / n) * n / dpi
    return (correct_roundoff(wnew, dpi, n), correct_roundoff(hnew, dpi, n))


# A registry for available MovieWriter classes
class MovieWriterRegistry(object):
    '''Registry of available writer classes by human readable name.'''
    def __init__(self):
        self.avail = dict()
        self._registered = dict()
        self._dirty = False

    def set_dirty(self):
        """Sets a flag to re-setup the writers."""
        self._dirty = True

    def register(self, name):
        """Decorator for registering a class under a name.

        Example use::

            @registry.register(name)
            class Foo:
                pass
        """
        def wrapper(writerClass):
            self._registered[name] = writerClass
            if writerClass.isAvailable():
                self.avail[name] = writerClass
            return writerClass
        return wrapper

    def ensure_not_dirty(self):
        """If dirty, reasks the writers if they are available"""
        if self._dirty:
            self.reset_available_writers()

    def reset_available_writers(self):
        """Reset the available state of all registered writers"""
        self.avail = {}
        for name, writerClass in self._registered.items():
            if writerClass.isAvailable():
                self.avail[name] = writerClass
        self._dirty = False

    def list(self):
        '''Get a list of available MovieWriters.'''
        self.ensure_not_dirty()
        return list(self.avail)

    def is_available(self, name):
        '''Check if given writer is available by name.

        Parameters
        ----------
        name : str

        Returns
        -------
        available : bool
        '''
        self.ensure_not_dirty()
        return name in self.avail

    def __getitem__(self, name):
        self.ensure_not_dirty()
        if not self.avail:
            raise RuntimeError("No MovieWriters available!")
        try:
            return self.avail[name]
        except KeyError:
            raise RuntimeError(
                'Requested MovieWriter ({}) not available'.format(name))


writers = MovieWriterRegistry()


class AbstractMovieWriter(six.with_metaclass(abc.ABCMeta)):
    '''
    Abstract base class for writing movies. Fundamentally, what a MovieWriter
    does is provide is a way to grab frames by calling grab_frame().

    setup() is called to start the process and finish() is called afterwards.

    This class is set up to provide for writing movie frame data to a pipe.
    saving() is provided as a context manager to facilitate this process as::

        with moviewriter.saving(fig, outfile='myfile.mp4', dpi=100):
            # Iterate over frames
            moviewriter.grab_frame(**savefig_kwargs)

    The use of the context manager ensures that setup() and finish() are
    performed as necessary.

    An instance of a concrete subclass of this class can be given as the
    ``writer`` argument of `Animation.save()`.
    '''

    @abc.abstractmethod
    def setup(self, fig, outfile, dpi=None):
        '''
        Perform setup for writing the movie file.

        Parameters
        ----------
        fig: `matplotlib.figure.Figure` instance
            The figure object that contains the information for frames
        outfile: string
            The filename of the resulting movie file
        dpi: int, optional
            The DPI (or resolution) for the file.  This controls the size
            in pixels of the resulting movie file. Default is ``fig.dpi``.
        '''

    @abc.abstractmethod
    def grab_frame(self, **savefig_kwargs):
        '''
        Grab the image information from the figure and save as a movie frame.

        All keyword arguments in savefig_kwargs are passed on to the `savefig`
        command that saves the figure.
        '''

    @abc.abstractmethod
    def finish(self):
        '''Finish any processing for writing the movie.'''

    @contextlib.contextmanager
    def saving(self, fig, outfile, dpi, *args, **kwargs):
        '''
        Context manager to facilitate writing the movie file.

        ``*args, **kw`` are any parameters that should be passed to `setup`.
        '''
        # This particular sequence is what contextlib.contextmanager wants
        self.setup(fig, outfile, dpi, *args, **kwargs)
        try:
            yield self
        finally:
            self.finish()


class MovieWriter(AbstractMovieWriter):
    '''Base class for writing movies.

    This class is set up to provide for writing movie frame data to a pipe.
    See examples for how to use these classes.

    Attributes
    ----------
    frame_format : str
        The format used in writing frame data, defaults to 'rgba'
    fig : `~matplotlib.figure.Figure`
        The figure to capture data from.
        This must be provided by the sub-classes.

    '''

    def __init__(self, fps=5, codec=None, bitrate=None, extra_args=None,
                 metadata=None):
        '''MovieWriter

        Parameters
        ----------
        fps: int
            Framerate for movie.
        codec: string or None, optional
            The codec to use. If ``None`` (the default) the ``animation.codec``
            rcParam is used.
        bitrate: int or None, optional
            The bitrate for the saved movie file, which is one way to control
            the output file size and quality. The default value is ``None``,
            which uses the ``animation.bitrate`` rcParam.  A value of -1
            implies that the bitrate should be determined automatically by the
            underlying utility.
        extra_args: list of strings or None, optional
            A list of extra string arguments to be passed to the underlying
            movie utility. The default is ``None``, which passes the additional
            arguments in the ``animation.extra_args`` rcParam.
        metadata: Dict[str, str] or None
            A dictionary of keys and values for metadata to include in the
            output file. Some keys that may be of use include:
            title, artist, genre, subject, copyright, srcform, comment.
        '''
        self.fps = fps
        self.frame_format = 'rgba'

        if codec is None:
            self.codec = rcParams['animation.codec']
        else:
            self.codec = codec

        if bitrate is None:
            self.bitrate = rcParams['animation.bitrate']
        else:
            self.bitrate = bitrate

        if extra_args is None:
            self.extra_args = list(rcParams[self.args_key])
        else:
            self.extra_args = extra_args

        if metadata is None:
            self.metadata = dict()
        else:
            self.metadata = metadata

    @property
    def frame_size(self):
        '''A tuple ``(width, height)`` in pixels of a movie frame.'''
        w, h = self.fig.get_size_inches()
        return int(w * self.dpi), int(h * self.dpi)

    def _adjust_frame_size(self):
        if self.codec == 'h264':
            wo, ho = self.fig.get_size_inches()
            w, h = adjusted_figsize(wo, ho, self.dpi, 2)
            if not (wo, ho) == (w, h):
                self.fig.set_size_inches(w, h, forward=True)
                _log.info('figure size (inches) has been adjusted '
                          'from %s x %s to %s x %s', wo, ho, w, h)
        else:
            w, h = self.fig.get_size_inches()
        _log.debug('frame size in pixels is %s x %s', *self.frame_size)
        return w, h

    def setup(self, fig, outfile, dpi=None):
        '''
        Perform setup for writing the movie file.

        Parameters
        ----------
        fig : matplotlib.figure.Figure
            The figure object that contains the information for frames
        outfile : string
            The filename of the resulting movie file
        dpi : int, optional
            The DPI (or resolution) for the file.  This controls the size
            in pixels of the resulting movie file. Default is fig.dpi.
        '''
        self.outfile = outfile
        self.fig = fig
        if dpi is None:
            dpi = self.fig.dpi
        self.dpi = dpi
        self._w, self._h = self._adjust_frame_size()

        # Run here so that grab_frame() can write the data to a pipe. This
        # eliminates the need for temp files.
        self._run()

    def _run(self):
        # Uses subprocess to call the program for assembling frames into a
        # movie file.  *args* returns the sequence of command line arguments
        # from a few configuration options.
        command = self._args()
        output = subprocess.PIPE
        _log.info('MovieWriter.run: running command: %s', command)
        self._proc = subprocess.Popen(command, shell=False,
                                      stdout=output, stderr=output,
                                      stdin=subprocess.PIPE,
                                      creationflags=subprocess_creation_flags)

    def finish(self):
        '''Finish any processing for writing the movie.'''
        self.cleanup()

    def grab_frame(self, **savefig_kwargs):
        '''
        Grab the image information from the figure and save as a movie frame.

        All keyword arguments in savefig_kwargs are passed on to the `savefig`
        command that saves the figure.
        '''
        _log.debug('MovieWriter.grab_frame: Grabbing frame.')
        try:
            # re-adjust the figure size in case it has been changed by the
            # user.  We must ensure that every frame is the same size or
            # the movie will not save correctly.
            self.fig.set_size_inches(self._w, self._h)
            # Tell the figure to save its data to the sink, using the
            # frame format and dpi.
            self.fig.savefig(self._frame_sink(), format=self.frame_format,
                             dpi=self.dpi, **savefig_kwargs)
        except (RuntimeError, IOError) as e:
            out, err = self._proc.communicate()
            _log.info('MovieWriter -- Error '
                           'running proc:\n%s\n%s' % (out, err))
            raise IOError('Error saving animation to file (cause: {0}) '
                          'Stdout: {1} StdError: {2}. It may help to re-run '
                          'with logging level set to '
                          'DEBUG.'.format(e, out, err))

    def _frame_sink(self):
        '''Returns the place to which frames should be written.'''
        return self._proc.stdin

    def _args(self):
        '''Assemble list of utility-specific command-line arguments.'''
        return NotImplementedError("args needs to be implemented by subclass.")

    def cleanup(self):
        '''Clean-up and collect the process used to write the movie file.'''
        out, err = self._proc.communicate()
        self._frame_sink().close()
        _log.debug('MovieWriter -- Command stdout:\n%s', out)
        _log.debug('MovieWriter -- Command stderr:\n%s', err)

    @classmethod
    def bin_path(cls):
        '''
        Returns the binary path to the commandline tool used by a specific
        subclass. This is a class method so that the tool can be looked for
        before making a particular MovieWriter subclass available.
        '''
        return str(rcParams[cls.exec_key])

    @classmethod
    def isAvailable(cls):
        '''
        Check to see if a MovieWriter subclass is actually available by
        running the commandline tool.
        '''
        bin_path = cls.bin_path()
        if not bin_path:
            return False
        try:
            p = subprocess.Popen(
                bin_path,
                shell=False,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                creationflags=subprocess_creation_flags)
            return cls._handle_subprocess(p)
        except OSError:
            return False

    @classmethod
    def _handle_subprocess(cls, process):
        process.communicate()
        return True


class FileMovieWriter(MovieWriter):
    '''`MovieWriter` for writing to individual files and stitching at the end.

    This must be sub-classed to be useful.
    '''
    def __init__(self, *args, **kwargs):
        MovieWriter.__init__(self, *args, **kwargs)
        self.frame_format = rcParams['animation.frame_format']

    def setup(self, fig, outfile, dpi=None, frame_prefix='_tmp',
              clear_temp=True):
        '''Perform setup for writing the movie file.

        Parameters
        ----------
        fig : matplotlib.figure.Figure
            The figure to grab the rendered frames from.
        outfile : str
            The filename of the resulting movie file.
        dpi : number, optional
            The dpi of the output file. This, with the figure size,
            controls the size in pixels of the resulting movie file.
            Default is fig.dpi.
        frame_prefix : str, optional
            The filename prefix to use for temporary files.  Defaults to
            ``'_tmp'``.
        clear_temp : bool, optional
            If the temporary files should be deleted after stitching
            the final result.  Setting this to ``False`` can be useful for
            debugging.  Defaults to ``True``.

        '''
        self.fig = fig
        self.outfile = outfile
        if dpi is None:
            dpi = self.fig.dpi
        self.dpi = dpi
        self._adjust_frame_size()

        self.clear_temp = clear_temp
        self.temp_prefix = frame_prefix
        self._frame_counter = 0  # used for generating sequential file names
        self._temp_names = list()
        self.fname_format_str = '%s%%07d.%s'

    @property
    def frame_format(self):
        '''
        Format (png, jpeg, etc.) to use for saving the frames, which can be
        decided by the individual subclasses.
        '''
        return self._frame_format

    @frame_format.setter
    def frame_format(self, frame_format):
        if frame_format in self.supported_formats:
            self._frame_format = frame_format
        else:
            self._frame_format = self.supported_formats[0]

    def _base_temp_name(self):
        # Generates a template name (without number) given the frame format
        # for extension and the prefix.
        return self.fname_format_str % (self.temp_prefix, self.frame_format)

    def _frame_sink(self):
        # Creates a filename for saving using the basename and the current
        # counter.
        fname = self._base_temp_name() % self._frame_counter

        # Save the filename so we can delete it later if necessary
        self._temp_names.append(fname)
        _log.debug('FileMovieWriter.frame_sink: saving frame %d to fname=%s',
                   self._frame_counter, fname)
        self._frame_counter += 1  # Ensures each created name is 'unique'

        # This file returned here will be closed once it's used by savefig()
        # because it will no longer be referenced and will be gc-ed.
        return open(fname, 'wb')

    def grab_frame(self, **savefig_kwargs):
        '''
        Grab the image information from the figure and save as a movie frame.
        All keyword arguments in savefig_kwargs are passed on to the `savefig`
        command that saves the figure.
        '''
        # Overloaded to explicitly close temp file.
        _log.debug('MovieWriter.grab_frame: Grabbing frame.')
        try:
            # Tell the figure to save its data to the sink, using the
            # frame format and dpi.
            with self._frame_sink() as myframesink:
                self.fig.savefig(myframesink, format=self.frame_format,
                                 dpi=self.dpi, **savefig_kwargs)

        except RuntimeError:
            out, err = self._proc.communicate()
            _log.info('MovieWriter -- Error '
                           'running proc:\n%s\n%s' % (out, err))
            raise

    def finish(self):
        # Call run here now that all frame grabbing is done. All temp files
        # are available to be assembled.
        self._run()
        MovieWriter.finish(self)  # Will call clean-up

        # Check error code for creating file here, since we just run
        # the process here, rather than having an open pipe.
        if self._proc.returncode:
            try:
                stdout = [s.decode() for s in self._proc._stdout_buff]
                stderr = [s.decode() for s in self._proc._stderr_buff]
                _log.info("MovieWriter.finish: stdout: %s", stdout)
                _log.info("MovieWriter.finish: stderr: %s", stderr)
            except Exception as e:
                pass
            raise RuntimeError('Error creating movie, return code: {}'
                               .format(self._proc.returncode))

    def cleanup(self):
        MovieWriter.cleanup(self)

        # Delete temporary files
        if self.clear_temp:
            _log.debug('MovieWriter: clearing temporary fnames=%s',
                       self._temp_names)
            for fname in self._temp_names:
                os.remove(fname)


@writers.register('pillow')
class PillowWriter(MovieWriter):
    @classmethod
    def isAvailable(cls):
        try:
            import PIL
        except ImportError:
            return False
        return True

    def __init__(self, *args, **kwargs):
        if kwargs.get("extra_args") is None:
            kwargs["extra_args"] = ()
        super(PillowWriter, self).__init__(*args, **kwargs)

    def setup(self, fig, outfile, dpi=None):
        self._frames = []
        self._outfile = outfile
        self._dpi = dpi
        self._fig = fig

    def grab_frame(self, **savefig_kwargs):
        from PIL import Image
        buf = BytesIO()
        self._fig.savefig(buf, **dict(savefig_kwargs, format="rgba"))
        renderer = self._fig.canvas.get_renderer()
        # Using frombuffer / getbuffer may be slightly more efficient, but
        # Py3-only.
        self._frames.append(Image.frombytes(
            "RGBA",
            (int(renderer.width), int(renderer.height)),
            buf.getvalue()))

    def finish(self):
        self._frames[0].save(
            self._outfile, save_all=True, append_images=self._frames[1:],
            duration=int(1000 / self.fps))


# Base class of ffmpeg information. Has the config keys and the common set
# of arguments that controls the *output* side of things.
class FFMpegBase(object):
    '''Mixin class for FFMpeg output.

    To be useful this must be multiply-inherited from with a
    `MovieWriterBase` sub-class.
    '''

    exec_key = 'animation.ffmpeg_path'
    args_key = 'animation.ffmpeg_args'

    @property
    def output_args(self):
        args = ['-vcodec', self.codec]
        # For h264, the default format is yuv444p, which is not compatible
        # with quicktime (and others). Specifying yuv420p fixes playback on
        # iOS,as well as HTML5 video in firefox and safari (on both Win and
        # OSX). Also fixes internet explorer. This is as of 2015/10/29.
        if self.codec == 'h264' and '-pix_fmt' not in self.extra_args:
            args.extend(['-pix_fmt', 'yuv420p'])
        # The %dk adds 'k' as a suffix so that ffmpeg treats our bitrate as in
        # kbps
        if self.bitrate > 0:
            args.extend(['-b', '%dk' % self.bitrate])
        if self.extra_args:
            args.extend(self.extra_args)
        for k, v in six.iteritems(self.metadata):
            args.extend(['-metadata', '%s=%s' % (k, v)])

        return args + ['-y', self.outfile]

    @classmethod
    def _handle_subprocess(cls, process):
        _, err = process.communicate()
        # Ubuntu 12.04 ships a broken ffmpeg binary which we shouldn't use
        # NOTE : when removed, remove the same method in AVConvBase.
        if 'Libav' in err.decode():
            return False
        return True


# Combine FFMpeg options with pipe-based writing
@writers.register('ffmpeg')
class FFMpegWriter(FFMpegBase, MovieWriter):
    '''Pipe-based ffmpeg writer.

    Frames are streamed directly to ffmpeg via a pipe and written in a single
    pass.
    '''
    def _args(self):
        # Returns the command line parameters for subprocess to use
        # ffmpeg to create a movie using a pipe.
        args = [self.bin_path(), '-f', 'rawvideo', '-vcodec', 'rawvideo',
                '-s', '%dx%d' % self.frame_size, '-pix_fmt', self.frame_format,
                '-r', str(self.fps)]
        # Logging is quieted because subprocess.PIPE has limited buffer size.
        # If you have a lot of frames in your animation and set logging to
        # DEBUG, you will have a buffer overrun.
        if (_log.getEffectiveLevel() > logging.DEBUG):
            args += ['-loglevel', 'quiet']
        args += ['-i', 'pipe:'] + self.output_args
        return args


# Combine FFMpeg options with temp file-based writing
@writers.register('ffmpeg_file')
class FFMpegFileWriter(FFMpegBase, FileMovieWriter):
    '''File-based ffmpeg writer.

    Frames are written to temporary files on disk and then stitched
    together at the end.

    '''
    supported_formats = ['png', 'jpeg', 'ppm', 'tiff', 'sgi', 'bmp',
                         'pbm', 'raw', 'rgba']

    def _args(self):
        # Returns the command line parameters for subprocess to use
        # ffmpeg to create a movie using a collection of temp images
        return [self.bin_path(), '-r', str(self.fps),
                '-i', self._base_temp_name(),
                '-vframes', str(self._frame_counter)] + self.output_args


# Base class of avconv information.  AVConv has identical arguments to
# FFMpeg
class AVConvBase(FFMpegBase):
    '''Mixin class for avconv output.

    To be useful this must be multiply-inherited from with a
    `MovieWriterBase` sub-class.
    '''

    exec_key = 'animation.avconv_path'
    args_key = 'animation.avconv_args'

    # NOTE : should be removed when the same method is removed in FFMpegBase.
    @classmethod
    def _handle_subprocess(cls, process):
        return MovieWriter._handle_subprocess(process)


# Combine AVConv options with pipe-based writing
@writers.register('avconv')
class AVConvWriter(AVConvBase, FFMpegWriter):
    '''Pipe-based avconv writer.

    Frames are streamed directly to avconv via a pipe and written in a single
    pass.
    '''


# Combine AVConv options with file-based writing
@writers.register('avconv_file')
class AVConvFileWriter(AVConvBase, FFMpegFileWriter):
    '''File-based avconv writer.

    Frames are written to temporary files on disk and then stitched
    together at the end.
    '''


# Base class for animated GIFs with convert utility
class ImageMagickBase(object):
    '''Mixin class for ImageMagick output.

    To be useful this must be multiply-inherited from with a
    `MovieWriterBase` sub-class.
    '''

    exec_key = 'animation.convert_path'
    args_key = 'animation.convert_args'

    @property
    def delay(self):
        return 100. / self.fps

    @property
    def output_args(self):
        return [self.outfile]

    @classmethod
    def _init_from_registry(cls):
        if sys.platform != 'win32' or rcParams[cls.exec_key] != 'convert':
            return
        from six.moves import winreg
        for flag in (0, winreg.KEY_WOW64_32KEY, winreg.KEY_WOW64_64KEY):
            try:
                hkey = winreg.OpenKeyEx(winreg.HKEY_LOCAL_MACHINE,
                                        'Software\\Imagemagick\\Current',
                                        0, winreg.KEY_QUERY_VALUE | flag)
                binpath = winreg.QueryValueEx(hkey, 'BinPath')[0]
                winreg.CloseKey(hkey)
                binpath += '\\convert.exe'
                break
            except Exception:
                binpath = ''
        rcParams[cls.exec_key] = rcParamsDefault[cls.exec_key] = binpath

    @classmethod
    def isAvailable(cls):
        '''
        Check to see if a ImageMagickWriter is actually available.

        Done by first checking the windows registry (if applicable) and then
        running the commandline tool.
        '''
        bin_path = cls.bin_path()
        if bin_path == "convert":
            cls._init_from_registry()
        return super(ImageMagickBase, cls).isAvailable()

ImageMagickBase._init_from_registry()


# Note: the base classes need to be in that order to get
# isAvailable() from ImageMagickBase called and not the
# one from MovieWriter. The latter is then called by the
# former.
@writers.register('imagemagick')
class ImageMagickWriter(ImageMagickBase, MovieWriter):
    '''Pipe-based animated gif.

    Frames are streamed directly to ImageMagick via a pipe and written
    in a single pass.

    '''
    def _args(self):
        return ([self.bin_path(),
                 '-size', '%ix%i' % self.frame_size, '-depth', '8',
                 '-delay', str(self.delay), '-loop', '0',
                 '%s:-' % self.frame_format]
                + self.output_args)


# Note: the base classes need to be in that order to get
# isAvailable() from ImageMagickBase called and not the
# one from MovieWriter. The latter is then called by the
# former.
@writers.register('imagemagick_file')
class ImageMagickFileWriter(ImageMagickBase, FileMovieWriter):
    '''File-based animated gif writer.

    Frames are written to temporary files on disk and then stitched
    together at the end.

    '''

    supported_formats = ['png', 'jpeg', 'ppm', 'tiff', 'sgi', 'bmp',
                         'pbm', 'raw', 'rgba']

    def _args(self):
        return ([self.bin_path(), '-delay', str(self.delay), '-loop', '0',
                 '%s*.%s' % (self.temp_prefix, self.frame_format)]
                + self.output_args)


# Taken directly from jakevdp's JSAnimation package at
# http://github.com/jakevdp/JSAnimation
def _included_frames(frame_list, frame_format):
    """frame_list should be a list of filenames"""
    return INCLUDED_FRAMES.format(Nframes=len(frame_list),
                                  frame_dir=os.path.dirname(frame_list[0]),
                                  frame_format=frame_format)


def _embedded_frames(frame_list, frame_format):
    """frame_list should be a list of base64-encoded png files"""
    template = '  frames[{0}] = "data:image/{1};base64,{2}"\n'
    return "\n" + "".join(
        template.format(i, frame_format, frame_data.replace('\n', '\\\n'))
        for i, frame_data in enumerate(frame_list))


@writers.register('html')
class HTMLWriter(FileMovieWriter):
    supported_formats = ['png', 'jpeg', 'tiff', 'svg']
    args_key = 'animation.html_args'

    @classmethod
    def isAvailable(cls):
        return True

    def __init__(self, fps=30, codec=None, bitrate=None, extra_args=None,
                 metadata=None, embed_frames=False, default_mode='loop',
                 embed_limit=None):
        self.embed_frames = embed_frames
        self.default_mode = default_mode.lower()

        # Save embed limit, which is given in MB
        if embed_limit is None:
            self._bytes_limit = rcParams['animation.embed_limit']
        else:
            self._bytes_limit = embed_limit

        # Convert from MB to bytes
        self._bytes_limit *= 1024 * 1024

        if self.default_mode not in ['loop', 'once', 'reflect']:
            self.default_mode = 'loop'
            _log.warning("unrecognized default_mode: using 'loop'")

        self._saved_frames = []
        self._total_bytes = 0
        self._hit_limit = False
        super(HTMLWriter, self).__init__(fps, codec, bitrate,
                                         extra_args, metadata)

    def setup(self, fig, outfile, dpi, frame_dir=None):
        root, ext = os.path.splitext(outfile)
        if ext not in ['.html', '.htm']:
            raise ValueError("outfile must be *.htm or *.html")

        if not self.embed_frames:
            if frame_dir is None:
                frame_dir = root + '_frames'
            if not os.path.exists(frame_dir):
                os.makedirs(frame_dir)
            frame_prefix = os.path.join(frame_dir, 'frame')
        else:
            frame_prefix = None

        super(HTMLWriter, self).setup(fig, outfile, dpi,
                                      frame_prefix, clear_temp=False)

    def grab_frame(self, **savefig_kwargs):
        if self.embed_frames:
            # Just stop processing if we hit the limit
            if self._hit_limit:
                return
            suffix = '.' + self.frame_format
            f = BytesIO()
            self.fig.savefig(f, format=self.frame_format,
                             dpi=self.dpi, **savefig_kwargs)
            imgdata64 = encodebytes(f.getvalue()).decode('ascii')
            self._total_bytes += len(imgdata64)
            if self._total_bytes >= self._bytes_limit:
                _log.warning(
                    "Animation size has reached %s bytes, exceeding the limit "
                    "of %s. If you're sure you want a larger animation "
                    "embedded, set the animation.embed_limit rc parameter to "
                    "a larger value (in MB). This and further frames will be "
                    "dropped.", self._total_bytes, self._bytes_limit)
                self._hit_limit = True
            else:
                self._saved_frames.append(imgdata64)
        else:
            return super(HTMLWriter, self).grab_frame(**savefig_kwargs)

    def _run(self):
        # make a duck-typed subprocess stand in
        # this is called by the MovieWriter base class, but not used here.
        class ProcessStandin(object):
            returncode = 0

            def communicate(self):
                return '', ''

        self._proc = ProcessStandin()

        # save the frames to an html file
        if self.embed_frames:
            fill_frames = _embedded_frames(self._saved_frames,
                                           self.frame_format)
        else:
            # temp names is filled by FileMovieWriter
            fill_frames = _included_frames(self._temp_names,
                                           self.frame_format)

        mode_dict = dict(once_checked='',
                         loop_checked='',
                         reflect_checked='')
        mode_dict[self.default_mode + '_checked'] = 'checked'

        interval = 1000 // self.fps

        with open(self.outfile, 'w') as of:
            of.write(JS_INCLUDE)
            of.write(DISPLAY_TEMPLATE.format(id=uuid.uuid4().hex,
                                             Nframes=len(self._temp_names),
                                             fill_frames=fill_frames,
                                             interval=interval,
                                             **mode_dict))


class Animation(object):
    '''This class wraps the creation of an animation using matplotlib.

    It is only a base class which should be subclassed to provide
    needed behavior.

    This class is not typically used directly.

    Parameters
    ----------
    fig : matplotlib.figure.Figure
       The figure object that is used to get draw, resize, and any
       other needed events.

    event_source : object, optional
       A class that can run a callback when desired events
       are generated, as well as be stopped and started.

       Examples include timers (see :class:`TimedAnimation`) and file
       system notifications.

    blit : bool, optional
       controls whether blitting is used to optimize drawing.  Defaults
       to ``False``.

    See Also
    --------
    FuncAnimation,  ArtistAnimation

    '''
    def __init__(self, fig, event_source=None, blit=False):
        self._fig = fig
        # Disables blitting for backends that don't support it.  This
        # allows users to request it if available, but still have a
        # fallback that works if it is not.
        self._blit = blit and fig.canvas.supports_blit

        # These are the basics of the animation.  The frame sequence represents
        # information for each frame of the animation and depends on how the
        # drawing is handled by the subclasses. The event source fires events
        # that cause the frame sequence to be iterated.
        self.frame_seq = self.new_frame_seq()
        self.event_source = event_source

        # Instead of starting the event source now, we connect to the figure's
        # draw_event, so that we only start once the figure has been drawn.
        self._first_draw_id = fig.canvas.mpl_connect('draw_event', self._start)

        # Connect to the figure's close_event so that we don't continue to
        # fire events and try to draw to a deleted figure.
        self._close_id = self._fig.canvas.mpl_connect('close_event',
                                                      self._stop)
        if self._blit:
            self._setup_blit()

    def _start(self, *args):
        '''
        Starts interactive animation. Adds the draw frame command to the GUI
        handler, calls show to start the event loop.
        '''
        # First disconnect our draw event handler
        self._fig.canvas.mpl_disconnect(self._first_draw_id)
        self._first_draw_id = None  # So we can check on save

        # Now do any initial draw
        self._init_draw()

        # Add our callback for stepping the animation and
        # actually start the event_source.
        self.event_source.add_callback(self._step)
        self.event_source.start()

    def _stop(self, *args):
        # On stop we disconnect all of our events.
        if self._blit:
            self._fig.canvas.mpl_disconnect(self._resize_id)
        self._fig.canvas.mpl_disconnect(self._close_id)
        self.event_source.remove_callback(self._step)
        self.event_source = None

    def save(self, filename, writer=None, fps=None, dpi=None, codec=None,
             bitrate=None, extra_args=None, metadata=None, extra_anim=None,
             savefig_kwargs=None):
        '''Saves a movie file by drawing every frame.

        Parameters
        ----------

        filename : str
            The output filename, e.g., :file:`mymovie.mp4`.

        writer : :class:`MovieWriter` or str, optional
            A `MovieWriter` instance to use or a key that identifies a
            class to use, such as 'ffmpeg'. If ``None``, defaults to
            :rc:`animation.writer`.

        fps : number, optional
           Frames per second in the movie. Defaults to ``None``, which will use
           the animation's specified interval to set the frames per second.

        dpi : number, optional
           Controls the dots per inch for the movie frames.  This combined with
           the figure's size in inches controls the size of the movie.  If
           ``None``, defaults to :rc:`savefig.dpi`.

        codec : str, optional
           The video codec to be used. Not all codecs are supported
           by a given :class:`MovieWriter`. If ``None``, default to
           :rc:`animation.codec`.

        bitrate : number, optional
           Specifies the number of bits used per second in the compressed
           movie, in kilobits per second. A higher number means a higher
           quality movie, but at the cost of increased file size. If ``None``,
           defaults to :rc:`animation.bitrate`.

        extra_args : list, optional
           List of extra string arguments to be passed to the underlying movie
           utility. If ``None``, defaults to :rc:`animation.extra_args`.

        metadata : Dict[str, str], optional
           Dictionary of keys and values for metadata to include in
           the output file. Some keys that may be of use include:
           title, artist, genre, subject, copyright, srcform, comment.

        extra_anim : list, optional
           Additional `Animation` objects that should be included
           in the saved movie file. These need to be from the same
           `matplotlib.figure.Figure` instance. Also, animation frames will
           just be simply combined, so there should be a 1:1 correspondence
           between the frames from the different animations.

        savefig_kwargs : dict, optional
           Is a dictionary containing keyword arguments to be passed
           on to the `savefig` command which is called repeatedly to
           save the individual frames.

        Notes
        -----
        fps, codec, bitrate, extra_args, metadata are used to
        construct a :class:`MovieWriter` instance and can only be
        passed if `writer` is a string.  If they are passed as
        non-`None` and ``writer`` is a :class:`MovieWriter`, a
        `RuntimeError` will be raised.

        '''
        # If the writer is None, use the rc param to find the name of the one
        # to use
        if writer is None:
            writer = rcParams['animation.writer']
        elif (not isinstance(writer, six.string_types) and
                any(arg is not None
                    for arg in (fps, codec, bitrate, extra_args, metadata))):
            raise RuntimeError('Passing in values for arguments '
                               'fps, codec, bitrate, extra_args, or metadata '
                               'is not supported when writer is an existing '
                               'MovieWriter instance. These should instead be '
                               'passed as arguments when creating the '
                               'MovieWriter instance.')

        if savefig_kwargs is None:
            savefig_kwargs = {}

        # Need to disconnect the first draw callback, since we'll be doing
        # draws. Otherwise, we'll end up starting the animation.
        if self._first_draw_id is not None:
            self._fig.canvas.mpl_disconnect(self._first_draw_id)
            reconnect_first_draw = True
        else:
            reconnect_first_draw = False

        if fps is None and hasattr(self, '_interval'):
            # Convert interval in ms to frames per second
            fps = 1000. / self._interval

        # Re-use the savefig DPI for ours if none is given
        if dpi is None:
            dpi = rcParams['savefig.dpi']
        if dpi == 'figure':
            dpi = self._fig.dpi

        if codec is None:
            codec = rcParams['animation.codec']

        if bitrate is None:
            bitrate = rcParams['animation.bitrate']

        all_anim = [self]
        if extra_anim is not None:
            all_anim.extend(anim
                            for anim
                            in extra_anim if anim._fig is self._fig)

        # If we have the name of a writer, instantiate an instance of the
        # registered class.
        if isinstance(writer, six.string_types):
            if writer in writers.avail:
                writer = writers[writer](fps, codec, bitrate,
                                         extra_args=extra_args,
                                         metadata=metadata)
            else:
                _log.warning("MovieWriter %s unavailable.", writer)

                try:
                    writer = writers[writers.list()[0]](fps, codec, bitrate,
                                                        extra_args=extra_args,
                                                        metadata=metadata)
                except IndexError:
                    raise ValueError("Cannot save animation: no writers are "
                                     "available. Please install ffmpeg to "
                                     "save animations.")
        _log.info('Animation.save using %s', type(writer))

        if 'bbox_inches' in savefig_kwargs:
            _log.warning("Warning: discarding the 'bbox_inches' argument in "
                         "'savefig_kwargs' as it may cause frame size "
                         "to vary, which is inappropriate for animation.")
            savefig_kwargs.pop('bbox_inches')

        # Create a new sequence of frames for saved data. This is different
        # from new_frame_seq() to give the ability to save 'live' generated
        # frame information to be saved later.
        # TODO: Right now, after closing the figure, saving a movie won't work
        # since GUI widgets are gone. Either need to remove extra code to
        # allow for this non-existent use case or find a way to make it work.
        with rc_context():
            if rcParams['savefig.bbox'] == 'tight':
                _log.info("Disabling savefig.bbox = 'tight', as it may cause "
                          "frame size to vary, which is inappropriate for "
                          "animation.")
                rcParams['savefig.bbox'] = None
            with writer.saving(self._fig, filename, dpi):
                for anim in all_anim:
                    # Clear the initial frame
                    anim._init_draw()
                for data in zip(*[a.new_saved_frame_seq() for a in all_anim]):
                    for anim, d in zip(all_anim, data):
                        # TODO: See if turning off blit is really necessary
                        anim._draw_next_frame(d, blit=False)
                    writer.grab_frame(**savefig_kwargs)

        # Reconnect signal for first draw if necessary
        if reconnect_first_draw:
            self._first_draw_id = self._fig.canvas.mpl_connect('draw_event',
                                                               self._start)

    def _step(self, *args):
        '''
        Handler for getting events. By default, gets the next frame in the
        sequence and hands the data off to be drawn.
        '''
        # Returns True to indicate that the event source should continue to
        # call _step, until the frame sequence reaches the end of iteration,
        # at which point False will be returned.
        try:
            framedata = next(self.frame_seq)
            self._draw_next_frame(framedata, self._blit)
            return True
        except StopIteration:
            return False

    def new_frame_seq(self):
        '''Creates a new sequence of frame information.'''
        # Default implementation is just an iterator over self._framedata
        return iter(self._framedata)

    def new_saved_frame_seq(self):
        '''Creates a new sequence of saved/cached frame information.'''
        # Default is the same as the regular frame sequence
        return self.new_frame_seq()

    def _draw_next_frame(self, framedata, blit):
        # Breaks down the drawing of the next frame into steps of pre- and
        # post- draw, as well as the drawing of the frame itself.
        self._pre_draw(framedata, blit)
        self._draw_frame(framedata)
        self._post_draw(framedata, blit)

    def _init_draw(self):
        # Initial draw to clear the frame. Also used by the blitting code
        # when a clean base is required.
        pass

    def _pre_draw(self, framedata, blit):
        # Perform any cleaning or whatnot before the drawing of the frame.
        # This default implementation allows blit to clear the frame.
        if blit:
            self._blit_clear(self._drawn_artists, self._blit_cache)

    def _draw_frame(self, framedata):
        # Performs actual drawing of the frame.
        raise NotImplementedError('Needs to be implemented by subclasses to'
                                  ' actually make an animation.')

    def _post_draw(self, framedata, blit):
        # After the frame is rendered, this handles the actual flushing of
        # the draw, which can be a direct draw_idle() or make use of the
        # blitting.
        if blit and self._drawn_artists:
            self._blit_draw(self._drawn_artists, self._blit_cache)
        else:
            self._fig.canvas.draw_idle()

    # The rest of the code in this class is to facilitate easy blitting
    def _blit_draw(self, artists, bg_cache):
        # Handles blitted drawing, which renders only the artists given instead
        # of the entire figure.
        updated_ax = []
        for a in artists:
            # If we haven't cached the background for this axes object, do
            # so now. This might not always be reliable, but it's an attempt
            # to automate the process.
            if a.axes not in bg_cache:
                bg_cache[a.axes] = a.figure.canvas.copy_from_bbox(a.axes.bbox)
            a.axes.draw_artist(a)
            updated_ax.append(a.axes)

        # After rendering all the needed artists, blit each axes individually.
        for ax in set(updated_ax):
            ax.figure.canvas.blit(ax.bbox)

    def _blit_clear(self, artists, bg_cache):
        # Get a list of the axes that need clearing from the artists that
        # have been drawn. Grab the appropriate saved background from the
        # cache and restore.
        axes = set(a.axes for a in artists)
        for a in axes:
            if a in bg_cache:
                a.figure.canvas.restore_region(bg_cache[a])

    def _setup_blit(self):
        # Setting up the blit requires: a cache of the background for the
        # axes
        self._blit_cache = dict()
        self._drawn_artists = []
        self._resize_id = self._fig.canvas.mpl_connect('resize_event',
                                                       self._handle_resize)
        self._post_draw(None, self._blit)

    def _handle_resize(self, *args):
        # On resize, we need to disable the resize event handling so we don't
        # get too many events. Also stop the animation events, so that
        # we're paused. Reset the cache and re-init. Set up an event handler
        # to catch once the draw has actually taken place.
        self._fig.canvas.mpl_disconnect(self._resize_id)
        self.event_source.stop()
        self._blit_cache.clear()
        self._init_draw()
        self._resize_id = self._fig.canvas.mpl_connect('draw_event',
                                                       self._end_redraw)

    def _end_redraw(self, evt):
        # Now that the redraw has happened, do the post draw flushing and
        # blit handling. Then re-enable all of the original events.
        self._post_draw(None, False)
        self.event_source.start()
        self._fig.canvas.mpl_disconnect(self._resize_id)
        self._resize_id = self._fig.canvas.mpl_connect('resize_event',
                                                       self._handle_resize)

    def to_html5_video(self, embed_limit=None):
        '''Returns animation as an HTML5 video tag.

        This saves the animation as an h264 video, encoded in base64
        directly into the HTML5 video tag. This respects the rc parameters
        for the writer as well as the bitrate. This also makes use of the
        ``interval`` to control the speed, and uses the ``repeat``
        parameter to decide whether to loop.
        '''
        VIDEO_TAG = r'''<video {size} {options}>
  <source type="video/mp4" src="data:video/mp4;base64,{video}">
  Your browser does not support the video tag.
</video>'''
        # Cache the rendering of the video as HTML
        if not hasattr(self, '_base64_video'):
            # Save embed limit, which is given in MB
            if embed_limit is None:
                embed_limit = rcParams['animation.embed_limit']

            # Convert from MB to bytes
            embed_limit *= 1024 * 1024

            # First write the video to a tempfile. Set delete to False
            # so we can re-open to read binary data.
            with tempfile.NamedTemporaryFile(suffix='.m4v',
                                             delete=False) as f:
                # We create a writer manually so that we can get the
                # appropriate size for the tag
                Writer = writers[rcParams['animation.writer']]
                writer = Writer(codec='h264',
                                bitrate=rcParams['animation.bitrate'],
                                fps=1000. / self._interval)
                self.save(f.name, writer=writer)

            # Now open and base64 encode
            with open(f.name, 'rb') as video:
                vid64 = encodebytes(video.read())
                vid_len = len(vid64)
                if vid_len >= embed_limit:
                    _log.warning(
                        "Animation movie is %s bytes, exceeding the limit of "
                        "%s. If you're sure you want a large animation "
                        "embedded, set the animation.embed_limit rc parameter "
                        "to a larger value (in MB).", vid_len, embed_limit)
                else:
                    self._base64_video = vid64.decode('ascii')
                    self._video_size = 'width="{}" height="{}"'.format(
                            *writer.frame_size)

            # Now we can remove
            os.remove(f.name)

        # If we exceeded the size, this attribute won't exist
        if hasattr(self, '_base64_video'):
            # Default HTML5 options are to autoplay and display video controls
            options = ['controls', 'autoplay']

            # If we're set to repeat, make it loop
            if hasattr(self, 'repeat') and self.repeat:
                options.append('loop')

            return VIDEO_TAG.format(video=self._base64_video,
                                    size=self._video_size,
                                    options=' '.join(options))
        else:
            return 'Video too large to embed.'

    def to_jshtml(self, fps=None, embed_frames=True, default_mode=None):
        """Generate HTML representation of the animation"""
        if fps is None and hasattr(self, '_interval'):
            # Convert interval in ms to frames per second
            fps = 1000 / self._interval

        # If we're not given a default mode, choose one base on the value of
        # the repeat attribute
        if default_mode is None:
            default_mode = 'loop' if self.repeat else 'once'

        if hasattr(self, "_html_representation"):
            return self._html_representation
        else:
            # Can't open a second time while opened on windows. So we avoid
            # deleting when closed, and delete manually later.
            with tempfile.NamedTemporaryFile(suffix='.html',
                                             delete=False) as f:
                self.save(f.name, writer=HTMLWriter(fps=fps,
                                                    embed_frames=embed_frames,
                                                    default_mode=default_mode))
            # Re-open and get content
            with open(f.name) as fobj:
                html = fobj.read()

            # Now we can delete
            os.remove(f.name)

            self._html_representation = html
            return html

    def _repr_html_(self):
        '''IPython display hook for rendering.'''
        fmt = rcParams['animation.html']
        if fmt == 'html5':
            return self.to_html5_video()
        elif fmt == 'jshtml':
            return self.to_jshtml()


class TimedAnimation(Animation):
    ''':class:`Animation` subclass for time-based animation.

    A new frame is drawn every *interval* milliseconds.

    Parameters
    ----------
    fig : matplotlib.figure.Figure
       The figure object that is used to get draw, resize, and any
       other needed events.

    interval : number, optional
       Delay between frames in milliseconds.  Defaults to 200.

    repeat_delay : number, optional
        If the animation in repeated, adds a delay in milliseconds
        before repeating the animation.  Defaults to ``None``.

    repeat : bool, optional
        Controls whether the animation should repeat when the sequence
        of frames is completed.  Defaults to ``True``.

    blit : bool, optional
        Controls whether blitting is used to optimize drawing.  Defaults
        to ``False``.

    '''
    def __init__(self, fig, interval=200, repeat_delay=None, repeat=True,
                 event_source=None, *args, **kwargs):
        # Store the timing information
        self._interval = interval
        self._repeat_delay = repeat_delay
        self.repeat = repeat

        # If we're not given an event source, create a new timer. This permits
        # sharing timers between animation objects for syncing animations.
        if event_source is None:
            event_source = fig.canvas.new_timer()
            event_source.interval = self._interval

        Animation.__init__(self, fig, event_source=event_source,
                           *args, **kwargs)

    def _step(self, *args):
        '''
        Handler for getting events.
        '''
        # Extends the _step() method for the Animation class.  If
        # Animation._step signals that it reached the end and we want to
        # repeat, we refresh the frame sequence and return True. If
        # _repeat_delay is set, change the event_source's interval to our loop
        # delay and set the callback to one which will then set the interval
        # back.
        still_going = Animation._step(self, *args)
        if not still_going and self.repeat:
            self._init_draw()
            self.frame_seq = self.new_frame_seq()
            if self._repeat_delay:
                self.event_source.remove_callback(self._step)
                self.event_source.add_callback(self._loop_delay)
                self.event_source.interval = self._repeat_delay
                return True
            else:
                return Animation._step(self, *args)
        else:
            return still_going

    def _stop(self, *args):
        # If we stop in the middle of a loop delay (which is relatively likely
        # given the potential pause here, remove the loop_delay callback as
        # well.
        self.event_source.remove_callback(self._loop_delay)
        Animation._stop(self)

    def _loop_delay(self, *args):
        # Reset the interval and change callbacks after the delay.
        self.event_source.remove_callback(self._loop_delay)
        self.event_source.interval = self._interval
        self.event_source.add_callback(self._step)
        Animation._step(self)


class ArtistAnimation(TimedAnimation):
    '''Animation using a fixed set of `Artist` objects.

    Before creating an instance, all plotting should have taken place
    and the relevant artists saved.

    Parameters
    ----------
    fig : matplotlib.figure.Figure
       The figure object that is used to get draw, resize, and any
       other needed events.

    artists : list
        Each list entry a collection of artists that represent what
        needs to be enabled on each frame. These will be disabled for
        other frames.

    interval : number, optional
       Delay between frames in milliseconds.  Defaults to 200.

    repeat_delay : number, optional
        If the animation in repeated, adds a delay in milliseconds
        before repeating the animation.  Defaults to ``None``.

    repeat : bool, optional
        Controls whether the animation should repeat when the sequence
        of frames is completed. Defaults to ``True``.

    blit : bool, optional
        Controls whether blitting is used to optimize drawing.  Defaults
        to ``False``.

    '''
    def __init__(self, fig, artists, *args, **kwargs):
        # Internal list of artists drawn in the most recent frame.
        self._drawn_artists = []

        # Use the list of artists as the framedata, which will be iterated
        # over by the machinery.
        self._framedata = artists
        TimedAnimation.__init__(self, fig, *args, **kwargs)

    def _init_draw(self):
        # Make all the artists involved in *any* frame invisible
        figs = set()
        for f in self.new_frame_seq():
            for artist in f:
                artist.set_visible(False)
                artist.set_animated(self._blit)
                # Assemble a list of unique figures that need flushing
                if artist.get_figure() not in figs:
                    figs.add(artist.get_figure())

        # Flush the needed figures
        for fig in figs:
            fig.canvas.draw_idle()

    def _pre_draw(self, framedata, blit):
        '''
        Clears artists from the last frame.
        '''
        if blit:
            # Let blit handle clearing
            self._blit_clear(self._drawn_artists, self._blit_cache)
        else:
            # Otherwise, make all the artists from the previous frame invisible
            for artist in self._drawn_artists:
                artist.set_visible(False)

    def _draw_frame(self, artists):
        # Save the artists that were passed in as framedata for the other
        # steps (esp. blitting) to use.
        self._drawn_artists = artists

        # Make all the artists from the current frame visible
        for artist in artists:
            artist.set_visible(True)


class FuncAnimation(TimedAnimation):
    '''
    Makes an animation by repeatedly calling a function ``func``.

    Parameters
    ----------
    fig : matplotlib.figure.Figure
       The figure object that is used to get draw, resize, and any
       other needed events.

    func : callable
       The function to call at each frame.  The first argument will
       be the next value in ``frames``.   Any additional positional
       arguments can be supplied via the ``fargs`` parameter.

       The required signature is::

          def func(frame, *fargs) -> iterable_of_artists:

    frames : iterable, int, generator function, or None, optional
        Source of data to pass ``func`` and each frame of the animation

        If an iterable, then simply use the values provided.  If the
        iterable has a length, it will override the ``save_count`` kwarg.

        If an integer, then equivalent to passing ``range(frames)``

        If a generator function, then must have the signature::

           def gen_function() -> obj:

        If ``None``, then equivalent to passing ``itertools.count``.

        In all of these cases, the values in *frames* is simply passed through
        to the user-supplied *func* and thus can be of any type.

    init_func : callable, optional
       A function used to draw a clear frame. If not given, the
       results of drawing from the first item in the frames sequence
       will be used. This function will be called once before the
       first frame.

       If ``blit == True``, ``init_func`` must return an iterable of artists
       to be re-drawn.

       The required signature is::

          def init_func() -> iterable_of_artists:

    fargs : tuple or None, optional
       Additional arguments to pass to each call to *func*.

    save_count : int, optional
       The number of values from *frames* to cache.

    interval : number, optional
       Delay between frames in milliseconds.  Defaults to 200.

    repeat_delay : number, optional
       If the animation in repeated, adds a delay in milliseconds
       before repeating the animation.  Defaults to ``None``.

    repeat : bool, optional
       Controls whether the animation should repeat when the sequence
       of frames is completed.  Defaults to ``True``.

    blit : bool, optional
       Controls whether blitting is used to optimize drawing.  Defaults
       to ``False``.

    '''
    def __init__(self, fig, func, frames=None, init_func=None, fargs=None,
                 save_count=None, **kwargs):
        if fargs:
            self._args = fargs
        else:
            self._args = ()
        self._func = func

        # Amount of framedata to keep around for saving movies. This is only
        # used if we don't know how many frames there will be: in the case
        # of no generator or in the case of a callable.
        self.save_count = save_count
        # Set up a function that creates a new iterable when needed. If nothing
        # is passed in for frames, just use itertools.count, which will just
        # keep counting from 0. A callable passed in for frames is assumed to
        # be a generator. An iterable will be used as is, and anything else
        # will be treated as a number of frames.
        if frames is None:
            self._iter_gen = itertools.count
        elif callable(frames):
            self._iter_gen = frames
        elif cbook.iterable(frames):
            self._iter_gen = lambda: iter(frames)
            if hasattr(frames, '__len__'):
                self.save_count = len(frames)
        else:
            self._iter_gen = lambda: iter(xrange(frames))
            self.save_count = frames

        if self.save_count is None:
            # If we're passed in and using the default, set save_count to 100.
            self.save_count = 100
        else:
            # itertools.islice returns an error when passed a numpy int instead
            # of a native python int (http://bugs.python.org/issue30537).
            # As a workaround, convert save_count to a native python int.
            self.save_count = int(self.save_count)

        self._init_func = init_func

        # Needs to be initialized so the draw functions work without checking
        self._save_seq = []

        TimedAnimation.__init__(self, fig, **kwargs)

        # Need to reset the saved seq, since right now it will contain data
        # for a single frame from init, which is not what we want.
        self._save_seq = []

    def new_frame_seq(self):
        # Use the generating function to generate a new frame sequence
        return self._iter_gen()

    def new_saved_frame_seq(self):
        # Generate an iterator for the sequence of saved data. If there are
        # no saved frames, generate a new frame sequence and take the first
        # save_count entries in it.
        if self._save_seq:
            # While iterating we are going to update _save_seq
            # so make a copy to safely iterate over
            self._old_saved_seq = list(self._save_seq)
            return iter(self._old_saved_seq)
        else:
            if self.save_count is not None:
                return itertools.islice(self.new_frame_seq(), self.save_count)

            else:
                frame_seq = self.new_frame_seq()

                def gen():
                    try:
                        for _ in range(100):
                            yield next(frame_seq)
                    except StopIteration:
                        pass
                    else:
                        cbook.warn_deprecated(
                            "2.2", "FuncAnimation.save has truncated your "
                            "animation to 100 frames.  In the future, no such "
                            "truncation will occur; please pass 'save_count' "
                            "accordingly.")

                return gen()

    def _init_draw(self):
        # Initialize the drawing either using the given init_func or by
        # calling the draw function with the first item of the frame sequence.
        # For blitting, the init_func should return a sequence of modified
        # artists.
        if self._init_func is None:
            self._draw_frame(next(self.new_frame_seq()))

        else:
            self._drawn_artists = self._init_func()
            if self._blit:
                if self._drawn_artists is None:
                    raise RuntimeError('The init_func must return a '
                                       'sequence of Artist objects.')
                for a in self._drawn_artists:
                    a.set_animated(self._blit)
        self._save_seq = []

    def _draw_frame(self, framedata):
        # Save the data for potential saving of movies.
        self._save_seq.append(framedata)

        # Make sure to respect save_count (keep only the last save_count
        # around)
        self._save_seq = self._save_seq[-self.save_count:]

        # Call the func with framedata and args. If blitting is desired,
        # func needs to return a sequence of any artists that were modified.
        self._drawn_artists = self._func(framedata, *self._args)
        if self._blit:
            if self._drawn_artists is None:
                raise RuntimeError('The animation function must return a '
                                   'sequence of Artist objects.')
            for a in self._drawn_artists:
                a.set_animated(self._blit)
