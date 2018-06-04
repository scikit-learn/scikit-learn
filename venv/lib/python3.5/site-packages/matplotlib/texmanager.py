r"""
This module supports embedded TeX expressions in matplotlib via dvipng
and dvips for the raster and postscript backends.  The tex and
dvipng/dvips information is cached in ~/.matplotlib/tex.cache for reuse between
sessions

Requirements:

* latex
* \*Agg backends: dvipng>=1.6
* PS backend: psfrag, dvips, and Ghostscript>=8.60

Backends:

* \*Agg
* PS
* PDF

For raster output, you can get RGBA numpy arrays from TeX expressions
as follows::

  texmanager = TexManager()
  s = ('\TeX\ is Number '
       '$\displaystyle\sum_{n=1}^\infty\frac{-e^{i\pi}}{2^n}$!')
  Z = texmanager.get_rgba(s, fontsize=12, dpi=80, rgb=(1,0,0))

To enable tex rendering of all text in your matplotlib figure, set
text.usetex in your matplotlibrc file or include these two lines in
your script::

  from matplotlib import rc
  rc('text', usetex=True)

"""

from __future__ import absolute_import, division, print_function

import six

import copy
import glob
import os
import shutil
import sys
import warnings
import logging

from hashlib import md5

import distutils.version
import numpy as np
import matplotlib as mpl
from matplotlib import rcParams
from matplotlib._png import read_png
from matplotlib.cbook import mkdirs, Locked
from matplotlib.compat.subprocess import subprocess, Popen, PIPE, STDOUT
import matplotlib.dviread as dviread
import re

_log = logging.getLogger(__name__)


@mpl.cbook.deprecated("2.1")
def dvipng_hack_alpha():
    try:
        p = Popen([str('dvipng'), '-version'], stdin=PIPE, stdout=PIPE,
                  stderr=STDOUT, close_fds=(sys.platform != 'win32'))
        stdout, stderr = p.communicate()
    except OSError:
        _log.info('No dvipng was found')
        return False
    lines = stdout.decode(sys.getdefaultencoding()).split('\n')
    for line in lines:
        if line.startswith('dvipng '):
            version = line.split()[-1]
            _log.info('Found dvipng version %s', version)
            version = distutils.version.LooseVersion(version)
            return version < distutils.version.LooseVersion('1.6')
    _log.info('Unexpected response from dvipng -version')
    return False


class TexManager(object):
    """
    Convert strings to dvi files using TeX, caching the results to a directory.
    """

    cachedir = mpl.get_cachedir()
    if cachedir is not None:
        texcache = os.path.join(cachedir, 'tex.cache')
        mkdirs(texcache)
    else:
        # Should only happen in a restricted environment (such as Google App
        # Engine). Deal with this gracefully by not creating a cache directory.
        texcache = None

    # Caches.
    rgba_arrayd = {}
    grey_arrayd = {}
    postscriptd = property(mpl.cbook.deprecated("2.2")(lambda self: {}))
    pscnt = property(mpl.cbook.deprecated("2.2")(lambda self: 0))

    serif = ('cmr', '')
    sans_serif = ('cmss', '')
    monospace = ('cmtt', '')
    cursive = ('pzc', r'\usepackage{chancery}')
    font_family = 'serif'
    font_families = ('serif', 'sans-serif', 'cursive', 'monospace')

    font_info = {
        'new century schoolbook': ('pnc', r'\renewcommand{\rmdefault}{pnc}'),
        'bookman': ('pbk', r'\renewcommand{\rmdefault}{pbk}'),
        'times': ('ptm', r'\usepackage{mathptmx}'),
        'palatino': ('ppl', r'\usepackage{mathpazo}'),
        'zapf chancery': ('pzc', r'\usepackage{chancery}'),
        'cursive': ('pzc', r'\usepackage{chancery}'),
        'charter': ('pch', r'\usepackage{charter}'),
        'serif': ('cmr', ''),
        'sans-serif': ('cmss', ''),
        'helvetica': ('phv', r'\usepackage{helvet}'),
        'avant garde': ('pag', r'\usepackage{avant}'),
        'courier': ('pcr', r'\usepackage{courier}'),
        'monospace': ('cmtt', ''),
        'computer modern roman': ('cmr', ''),
        'computer modern sans serif': ('cmss', ''),
        'computer modern typewriter': ('cmtt', '')}

    _rc_cache = None
    _rc_cache_keys = (('text.latex.preamble', ) +
                      tuple(['font.' + n for n in ('family', ) +
                             font_families]))

    def __init__(self):

        if self.texcache is None:
            raise RuntimeError('Cannot create TexManager, as there is no '
                               'cache directory available')

        mkdirs(self.texcache)
        ff = rcParams['font.family']
        if len(ff) == 1 and ff[0].lower() in self.font_families:
            self.font_family = ff[0].lower()
        elif (isinstance(ff, six.string_types)
              and ff.lower() in self.font_families):
            self.font_family = ff.lower()
        else:
            _log.info('font.family must be one of (%s) when text.usetex is '
                      'True. serif will be used by default.',
                      ', '.join(self.font_families))
            self.font_family = 'serif'

        fontconfig = [self.font_family]
        for font_family in self.font_families:
            font_family_attr = font_family.replace('-', '_')
            for font in rcParams['font.' + font_family]:
                if font.lower() in self.font_info:
                    setattr(self, font_family_attr,
                            self.font_info[font.lower()])
                    _log.debug('family: %s, font: %s, info: %s',
                               font_family, font, self.font_info[font.lower()])
                    break
                else:
                    _log.debug('%s font is not compatible with usetex.',
                               font_family)
            else:
                _log.info('No LaTeX-compatible font found for the %s font '
                          'family in rcParams. Using default.', font_family)
                setattr(self, font_family_attr, self.font_info[font_family])
            fontconfig.append(getattr(self, font_family_attr)[0])
        # Add a hash of the latex preamble to self._fontconfig so that the
        # correct png is selected for strings rendered with same font and dpi
        # even if the latex preamble changes within the session
        preamble_bytes = self.get_custom_preamble().encode('utf-8')
        fontconfig.append(md5(preamble_bytes).hexdigest())
        self._fontconfig = ''.join(fontconfig)

        # The following packages and commands need to be included in the latex
        # file's preamble:
        cmd = [self.serif[1], self.sans_serif[1], self.monospace[1]]
        if self.font_family == 'cursive':
            cmd.append(self.cursive[1])
        self._font_preamble = '\n'.join(
            [r'\usepackage{type1cm}'] + cmd + [r'\usepackage{textcomp}'])

    def get_basefile(self, tex, fontsize, dpi=None):
        """
        Return a filename based on a hash of the string, fontsize, and dpi.
        """
        s = ''.join([tex, self.get_font_config(), '%f' % fontsize,
                     self.get_custom_preamble(), str(dpi or '')])
        return os.path.join(self.texcache, md5(s.encode('utf-8')).hexdigest())

    def get_font_config(self):
        """Reinitializes self if relevant rcParams on have changed."""
        if self._rc_cache is None:
            self._rc_cache = dict.fromkeys(self._rc_cache_keys)
        changed = [par for par in self._rc_cache_keys
                   if rcParams[par] != self._rc_cache[par]]
        if changed:
            _log.debug('following keys changed: %s', changed)
            for k in changed:
                _log.debug('%-20s: %-10s -> %-10s',
                           k, self._rc_cache[k], rcParams[k])
                # deepcopy may not be necessary, but feels more future-proof
                self._rc_cache[k] = copy.deepcopy(rcParams[k])
            _log.debug('RE-INIT\nold fontconfig: %s', self._fontconfig)
            self.__init__()
        _log.debug('fontconfig: %s', self._fontconfig)
        return self._fontconfig

    def get_font_preamble(self):
        """
        Return a string containing font configuration for the tex preamble.
        """
        return self._font_preamble

    def get_custom_preamble(self):
        """Return a string containing user additions to the tex preamble."""
        return '\n'.join(rcParams['text.latex.preamble'])

    def make_tex(self, tex, fontsize):
        """
        Generate a tex file to render the tex string at a specific font size.

        Return the file name.
        """
        basefile = self.get_basefile(tex, fontsize)
        texfile = '%s.tex' % basefile
        custom_preamble = self.get_custom_preamble()
        fontcmd = {'sans-serif': r'{\sffamily %s}',
                   'monospace': r'{\ttfamily %s}'}.get(self.font_family,
                                                       r'{\rmfamily %s}')
        tex = fontcmd % tex

        if rcParams['text.latex.unicode']:
            unicode_preamble = r"""
\usepackage{ucs}
\usepackage[utf8x]{inputenc}"""
        else:
            unicode_preamble = ''

        s = r"""
\documentclass{article}
%s
%s
%s
\usepackage[papersize={72in,72in},body={70in,70in},margin={1in,1in}]{geometry}
\pagestyle{empty}
\begin{document}
\fontsize{%f}{%f}%s
\end{document}
""" % (self._font_preamble, unicode_preamble, custom_preamble,
       fontsize, fontsize * 1.25, tex)
        with open(texfile, 'wb') as fh:
            if rcParams['text.latex.unicode']:
                fh.write(s.encode('utf8'))
            else:
                try:
                    fh.write(s.encode('ascii'))
                except UnicodeEncodeError as err:
                    _log.info("You are using unicode and latex, but have not "
                              "enabled the 'text.latex.unicode' rcParam.")
                    raise

        return texfile

    _re_vbox = re.compile(
        r"MatplotlibBox:\(([\d.]+)pt\+([\d.]+)pt\)x([\d.]+)pt")

    def make_tex_preview(self, tex, fontsize):
        """
        Generate a tex file to render the tex string at a specific font size.

        It uses the preview.sty to determine the dimension (width, height,
        descent) of the output.

        Return the file name.
        """
        basefile = self.get_basefile(tex, fontsize)
        texfile = '%s.tex' % basefile
        custom_preamble = self.get_custom_preamble()
        fontcmd = {'sans-serif': r'{\sffamily %s}',
                   'monospace': r'{\ttfamily %s}'}.get(self.font_family,
                                                       r'{\rmfamily %s}')
        tex = fontcmd % tex

        if rcParams['text.latex.unicode']:
            unicode_preamble = r"""
\usepackage{ucs}
\usepackage[utf8x]{inputenc}"""
        else:
            unicode_preamble = ''

        # newbox, setbox, immediate, etc. are used to find the box
        # extent of the rendered text.

        s = r"""
\documentclass{article}
%s
%s
%s
\usepackage[active,showbox,tightpage]{preview}
\usepackage[papersize={72in,72in},body={70in,70in},margin={1in,1in}]{geometry}

%% we override the default showbox as it is treated as an error and makes
%% the exit status not zero
\def\showbox#1%%
{\immediate\write16{MatplotlibBox:(\the\ht#1+\the\dp#1)x\the\wd#1}}

\begin{document}
\begin{preview}
{\fontsize{%f}{%f}%s}
\end{preview}
\end{document}
""" % (self._font_preamble, unicode_preamble, custom_preamble,
       fontsize, fontsize * 1.25, tex)
        with open(texfile, 'wb') as fh:
            if rcParams['text.latex.unicode']:
                fh.write(s.encode('utf8'))
            else:
                try:
                    fh.write(s.encode('ascii'))
                except UnicodeEncodeError as err:
                    _log.info("You are using unicode and latex, but have not "
                              "enabled the 'text.latex.unicode' rcParam.")
                    raise

        return texfile

    def _run_checked_subprocess(self, command, tex):
        _log.debug(command)
        try:
            report = subprocess.check_output(command,
                                             cwd=self.texcache,
                                             stderr=subprocess.STDOUT)
        except subprocess.CalledProcessError as exc:
            raise RuntimeError(
                '{prog} was not able to process the following string:\n'
                '{tex!r}\n\n'
                'Here is the full report generated by {prog}:\n'
                '{exc}\n\n'.format(
                    prog=command[0],
                    tex=tex.encode('unicode_escape'),
                    exc=exc.output.decode('utf-8')))
        _log.debug(report)
        return report

    def make_dvi(self, tex, fontsize):
        """
        Generate a dvi file containing latex's layout of tex string.

        Return the file name.
        """

        if rcParams['text.latex.preview']:
            return self.make_dvi_preview(tex, fontsize)

        basefile = self.get_basefile(tex, fontsize)
        dvifile = '%s.dvi' % basefile
        if not os.path.exists(dvifile):
            texfile = self.make_tex(tex, fontsize)
            with Locked(self.texcache):
                self._run_checked_subprocess(
                    ["latex", "-interaction=nonstopmode", "--halt-on-error",
                     texfile], tex)
            for fname in glob.glob(basefile + '*'):
                if not fname.endswith(('dvi', 'tex')):
                    try:
                        os.remove(fname)
                    except OSError:
                        pass

        return dvifile

    def make_dvi_preview(self, tex, fontsize):
        """
        Generate a dvi file containing latex's layout of tex string.

        It calls make_tex_preview() method and store the size information
        (width, height, descent) in a separate file.

        Return the file name.
        """
        basefile = self.get_basefile(tex, fontsize)
        dvifile = '%s.dvi' % basefile
        baselinefile = '%s.baseline' % basefile

        if not os.path.exists(dvifile) or not os.path.exists(baselinefile):
            texfile = self.make_tex_preview(tex, fontsize)
            report = self._run_checked_subprocess(
                ["latex", "-interaction=nonstopmode", "--halt-on-error",
                 texfile], tex)

            # find the box extent information in the latex output
            # file and store them in ".baseline" file
            m = TexManager._re_vbox.search(report.decode("utf-8"))
            with open(basefile + '.baseline', "w") as fh:
                fh.write(" ".join(m.groups()))

            for fname in glob.glob(basefile + '*'):
                if not fname.endswith(('dvi', 'tex', 'baseline')):
                    try:
                        os.remove(fname)
                    except OSError:
                        pass

        return dvifile

    def make_png(self, tex, fontsize, dpi):
        """
        Generate a png file containing latex's rendering of tex string.

        Return the file name.
        """
        basefile = self.get_basefile(tex, fontsize, dpi)
        pngfile = '%s.png' % basefile
        # see get_rgba for a discussion of the background
        if not os.path.exists(pngfile):
            dvifile = self.make_dvi(tex, fontsize)
            self._run_checked_subprocess(
                ["dvipng", "-bg", "Transparent", "-D", str(dpi),
                 "-T", "tight", "-o", pngfile, dvifile], tex)
        return pngfile

    @mpl.cbook.deprecated("2.2")
    def make_ps(self, tex, fontsize):
        """
        Generate a postscript file containing latex's rendering of tex string.

        Return the file name.
        """
        basefile = self.get_basefile(tex, fontsize)
        psfile = '%s.epsf' % basefile
        if not os.path.exists(psfile):
            dvifile = self.make_dvi(tex, fontsize)
            self._run_checked_subprocess(
                ["dvips", "-q", "-E", "-o", psfile, dvifile], tex)
        return psfile

    @mpl.cbook.deprecated("2.2")
    def get_ps_bbox(self, tex, fontsize):
        """
        Return a list of PS bboxes for latex's rendering of the tex string.
        """
        psfile = self.make_ps(tex, fontsize)
        with open(psfile) as ps:
            for line in ps:
                if line.startswith('%%BoundingBox:'):
                    return [int(val) for val in line.split()[1:]]
        raise RuntimeError('Could not parse %s' % psfile)

    def get_grey(self, tex, fontsize=None, dpi=None):
        """Return the alpha channel."""
        key = tex, self.get_font_config(), fontsize, dpi
        alpha = self.grey_arrayd.get(key)
        if alpha is None:
            pngfile = self.make_png(tex, fontsize, dpi)
            X = read_png(os.path.join(self.texcache, pngfile))
            self.grey_arrayd[key] = alpha = X[:, :, -1]
        return alpha

    def get_rgba(self, tex, fontsize=None, dpi=None, rgb=(0, 0, 0)):
        """Return latex's rendering of the tex string as an rgba array."""
        if not fontsize:
            fontsize = rcParams['font.size']
        if not dpi:
            dpi = rcParams['savefig.dpi']
        r, g, b = rgb
        key = tex, self.get_font_config(), fontsize, dpi, tuple(rgb)
        Z = self.rgba_arrayd.get(key)

        if Z is None:
            alpha = self.get_grey(tex, fontsize, dpi)
            Z = np.dstack([r, g, b, alpha])
            self.rgba_arrayd[key] = Z

        return Z

    def get_text_width_height_descent(self, tex, fontsize, renderer=None):
        """Return width, height and descent of the text."""
        if tex.strip() == '':
            return 0, 0, 0

        dpi_fraction = renderer.points_to_pixels(1.) if renderer else 1

        if rcParams['text.latex.preview']:
            # use preview.sty
            basefile = self.get_basefile(tex, fontsize)
            baselinefile = '%s.baseline' % basefile

            if not os.path.exists(baselinefile):
                dvifile = self.make_dvi_preview(tex, fontsize)

            with open(baselinefile) as fh:
                l = fh.read().split()
            height, depth, width = [float(l1) * dpi_fraction for l1 in l]
            return width, height + depth, depth

        else:
            # use dviread. It sometimes returns a wrong descent.
            dvifile = self.make_dvi(tex, fontsize)
            with dviread.Dvi(dvifile, 72 * dpi_fraction) as dvi:
                page = next(iter(dvi))
            # A total height (including the descent) needs to be returned.
            return page.width, page.height + page.descent, page.descent
