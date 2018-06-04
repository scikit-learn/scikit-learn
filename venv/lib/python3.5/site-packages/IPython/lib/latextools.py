# -*- coding: utf-8 -*-
"""Tools for handling LaTeX."""

# Copyright (c) IPython Development Team.
# Distributed under the terms of the Modified BSD License.

from io import BytesIO, open
import os
import tempfile
import shutil
import subprocess
from base64 import encodebytes

from IPython.utils.process import find_cmd, FindCmdError
from traitlets.config import get_config
from traitlets.config.configurable import SingletonConfigurable
from traitlets import List, Bool, Unicode
from IPython.utils.py3compat import cast_unicode


class LaTeXTool(SingletonConfigurable):
    """An object to store configuration of the LaTeX tool."""
    def _config_default(self):
        return get_config()
    
    backends = List(
        Unicode(), ["matplotlib", "dvipng"],
        help="Preferred backend to draw LaTeX math equations. "
        "Backends in the list are checked one by one and the first "
        "usable one is used.  Note that `matplotlib` backend "
        "is usable only for inline style equations.  To draw  "
        "display style equations, `dvipng` backend must be specified. ",
        # It is a List instead of Enum, to make configuration more
        # flexible.  For example, to use matplotlib mainly but dvipng
        # for display style, the default ["matplotlib", "dvipng"] can
        # be used.  To NOT use dvipng so that other repr such as
        # unicode pretty printing is used, you can use ["matplotlib"].
        ).tag(config=True)

    use_breqn = Bool(
        True,
        help="Use breqn.sty to automatically break long equations. "
        "This configuration takes effect only for dvipng backend.",
        ).tag(config=True)

    packages = List(
        ['amsmath', 'amsthm', 'amssymb', 'bm'],
        help="A list of packages to use for dvipng backend. "
        "'breqn' will be automatically appended when use_breqn=True.",
        ).tag(config=True)

    preamble = Unicode(
        help="Additional preamble to use when generating LaTeX source "
        "for dvipng backend.",
        ).tag(config=True)


def latex_to_png(s, encode=False, backend=None, wrap=False):
    """Render a LaTeX string to PNG.

    Parameters
    ----------
    s : str
        The raw string containing valid inline LaTeX.
    encode : bool, optional
        Should the PNG data base64 encoded to make it JSON'able.
    backend : {matplotlib, dvipng}
        Backend for producing PNG data.
    wrap : bool
        If true, Automatically wrap `s` as a LaTeX equation.

    None is returned when the backend cannot be used.

    """
    s = cast_unicode(s)
    allowed_backends = LaTeXTool.instance().backends
    if backend is None:
        backend = allowed_backends[0]
    if backend not in allowed_backends:
        return None
    if backend == 'matplotlib':
        f = latex_to_png_mpl
    elif backend == 'dvipng':
        f = latex_to_png_dvipng
    else:
        raise ValueError('No such backend {0}'.format(backend))
    bin_data = f(s, wrap)
    if encode and bin_data:
        bin_data = encodebytes(bin_data)
    return bin_data


def latex_to_png_mpl(s, wrap):
    try:
        from matplotlib import mathtext
        from pyparsing import ParseFatalException
    except ImportError:
        return None

    # mpl mathtext doesn't support display math, force inline
    s = s.replace('$$', '$')
    if wrap:
        s = u'${0}$'.format(s)

    try:
        mt = mathtext.MathTextParser('bitmap')
        f = BytesIO()
        mt.to_png(f, s, fontsize=12)
        return f.getvalue()
    except (ValueError, RuntimeError, ParseFatalException):
        return None


def latex_to_png_dvipng(s, wrap):
    try:
        find_cmd('latex')
        find_cmd('dvipng')
    except FindCmdError:
        return None
    try:
        workdir = tempfile.mkdtemp()
        tmpfile = os.path.join(workdir, "tmp.tex")
        dvifile = os.path.join(workdir, "tmp.dvi")
        outfile = os.path.join(workdir, "tmp.png")

        with open(tmpfile, "w", encoding='utf8') as f:
            f.writelines(genelatex(s, wrap))

        with open(os.devnull, 'wb') as devnull:
            subprocess.check_call(
                ["latex", "-halt-on-error", "-interaction", "batchmode", tmpfile],
                cwd=workdir, stdout=devnull, stderr=devnull)

            subprocess.check_call(
                ["dvipng", "-T", "tight", "-x", "1500", "-z", "9",
                 "-bg", "transparent", "-o", outfile, dvifile], cwd=workdir,
                stdout=devnull, stderr=devnull)

        with open(outfile, "rb") as f:
            return f.read()
    except subprocess.CalledProcessError:
        return None
    finally:
        shutil.rmtree(workdir)


def kpsewhich(filename):
    """Invoke kpsewhich command with an argument `filename`."""
    try:
        find_cmd("kpsewhich")
        proc = subprocess.Popen(
            ["kpsewhich", filename],
            stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        (stdout, stderr) = proc.communicate()
        return stdout.strip().decode('utf8', 'replace')
    except FindCmdError:
        pass


def genelatex(body, wrap):
    """Generate LaTeX document for dvipng backend."""
    lt = LaTeXTool.instance()
    breqn = wrap and lt.use_breqn and kpsewhich("breqn.sty")
    yield r'\documentclass{article}'
    packages = lt.packages
    if breqn:
        packages = packages + ['breqn']
    for pack in packages:
        yield r'\usepackage{{{0}}}'.format(pack)
    yield r'\pagestyle{empty}'
    if lt.preamble:
        yield lt.preamble
    yield r'\begin{document}'
    if breqn:
        yield r'\begin{dmath*}'
        yield body
        yield r'\end{dmath*}'
    elif wrap:
        yield u'$${0}$$'.format(body)
    else:
        yield body
    yield u'\end{document}'


_data_uri_template_png = u"""<img src="data:image/png;base64,%s" alt=%s />"""

def latex_to_html(s, alt='image'):
    """Render LaTeX to HTML with embedded PNG data using data URIs.

    Parameters
    ----------
    s : str
        The raw string containing valid inline LateX.
    alt : str
        The alt text to use for the HTML.
    """
    base64_data = latex_to_png(s, encode=True).decode('ascii')
    if base64_data:
        return _data_uri_template_png  % (base64_data, alt)


