"""
    sphinx.highlighting
    ~~~~~~~~~~~~~~~~~~~

    Highlight code blocks using Pygments.

    :copyright: Copyright 2007-2022 by the Sphinx team, see AUTHORS.
    :license: BSD, see LICENSE for details.
"""

from functools import partial
from importlib import import_module
from typing import Any, Dict

from packaging import version
from pygments import __version__ as pygmentsversion
from pygments import highlight
from pygments.filters import ErrorToken
from pygments.formatter import Formatter
from pygments.formatters import HtmlFormatter, LatexFormatter
from pygments.lexer import Lexer
from pygments.lexers import (CLexer, Python3Lexer, PythonConsoleLexer, PythonLexer, RstLexer,
                             TextLexer, get_lexer_by_name, guess_lexer)
from pygments.style import Style
from pygments.styles import get_style_by_name
from pygments.util import ClassNotFound

from sphinx.locale import __
from sphinx.pygments_styles import NoneStyle, SphinxStyle
from sphinx.util import logging, texescape

logger = logging.getLogger(__name__)

lexers: Dict[str, Lexer] = {}
lexer_classes: Dict[str, Lexer] = {
    'none': partial(TextLexer, stripnl=False),
    'python': partial(PythonLexer, stripnl=False),
    'python3': partial(Python3Lexer, stripnl=False),
    'pycon': partial(PythonConsoleLexer, stripnl=False),
    'pycon3': partial(PythonConsoleLexer, python3=True, stripnl=False),
    'rest': partial(RstLexer, stripnl=False),
    'c': partial(CLexer, stripnl=False),
}


escape_hl_chars = {ord('\\'): '\\PYGZbs{}',
                   ord('{'): '\\PYGZob{}',
                   ord('}'): '\\PYGZcb{}'}

# used if Pygments is available
# use textcomp quote to get a true single quote
_LATEX_ADD_STYLES = r'''
\renewcommand\PYGZsq{\textquotesingle}
'''
# fix extra space between lines when Pygments highlighting uses \fcolorbox
# add a {..} to limit \fboxsep scope, and force \fcolorbox use correct value
# cf pygments #1708 which makes this unneeded for Pygments > 2.7.4
_LATEX_ADD_STYLES_FIXPYG = r'''
\makeatletter
% fix for Pygments <= 2.7.4
\let\spx@original@fcolorbox\fcolorbox
\def\spx@fixpyg@fcolorbox{\fboxsep-\fboxrule\spx@original@fcolorbox}
\def\PYG#1#2{\PYG@reset\PYG@toks#1+\relax+%
             {\let\fcolorbox\spx@fixpyg@fcolorbox\PYG@do{#2}}}
\makeatother
'''
if version.parse(pygmentsversion).release <= (2, 7, 4):
    _LATEX_ADD_STYLES += _LATEX_ADD_STYLES_FIXPYG


class PygmentsBridge:
    # Set these attributes if you want to have different Pygments formatters
    # than the default ones.
    html_formatter = HtmlFormatter
    latex_formatter = LatexFormatter

    def __init__(self, dest: str = 'html', stylename: str = 'sphinx',
                 latex_engine: str = None) -> None:
        self.dest = dest
        self.latex_engine = latex_engine

        style = self.get_style(stylename)
        self.formatter_args: Dict[str, Any] = {'style': style}
        if dest == 'html':
            self.formatter = self.html_formatter
        else:
            self.formatter = self.latex_formatter
            self.formatter_args['commandprefix'] = 'PYG'

    def get_style(self, stylename: str) -> Style:
        if stylename is None or stylename == 'sphinx':
            return SphinxStyle
        elif stylename == 'none':
            return NoneStyle
        elif '.' in stylename:
            module, stylename = stylename.rsplit('.', 1)
            return getattr(import_module(module), stylename)
        else:
            return get_style_by_name(stylename)

    def get_formatter(self, **kwargs: Any) -> Formatter:
        kwargs.update(self.formatter_args)
        return self.formatter(**kwargs)

    def get_lexer(self, source: str, lang: str, opts: Dict = None,
                  force: bool = False, location: Any = None) -> Lexer:
        if not opts:
            opts = {}

        # find out which lexer to use
        if lang in ('py', 'python'):
            if source.startswith('>>>'):
                # interactive session
                lang = 'pycon'
            else:
                lang = 'python'
        elif lang in ('py3', 'python3', 'default'):
            if source.startswith('>>>'):
                lang = 'pycon3'
            else:
                lang = 'python3'

        if lang in lexers:
            # just return custom lexers here (without installing raiseonerror filter)
            return lexers[lang]
        elif lang in lexer_classes:
            lexer = lexer_classes[lang](**opts)
        else:
            try:
                if lang == 'guess':
                    lexer = guess_lexer(source, **opts)
                else:
                    lexer = get_lexer_by_name(lang, **opts)
            except ClassNotFound:
                logger.warning(__('Pygments lexer name %r is not known'), lang,
                               location=location)
                lexer = lexer_classes['none'](**opts)

        if not force:
            lexer.add_filter('raiseonerror')

        return lexer

    def highlight_block(self, source: str, lang: str, opts: Dict = None,
                        force: bool = False, location: Any = None, **kwargs: Any) -> str:
        if not isinstance(source, str):
            source = source.decode()

        lexer = self.get_lexer(source, lang, opts, force, location)

        # highlight via Pygments
        formatter = self.get_formatter(**kwargs)
        try:
            hlsource = highlight(source, lexer, formatter)
        except ErrorToken:
            # this is most probably not the selected language,
            # so let it pass unhighlighted
            if lang == 'default':
                pass  # automatic highlighting failed.
            else:
                logger.warning(__('Could not lex literal_block as "%s". '
                                  'Highlighting skipped.'), lang,
                               type='misc', subtype='highlighting_failure',
                               location=location)
            lexer = self.get_lexer(source, 'none', opts, force, location)
            hlsource = highlight(source, lexer, formatter)

        if self.dest == 'html':
            return hlsource
        else:
            # MEMO: this is done to escape Unicode chars with non-Unicode engines
            return texescape.hlescape(hlsource, self.latex_engine)

    def get_stylesheet(self) -> str:
        formatter = self.get_formatter()
        if self.dest == 'html':
            return formatter.get_style_defs('.highlight')
        else:
            return formatter.get_style_defs() + _LATEX_ADD_STYLES
