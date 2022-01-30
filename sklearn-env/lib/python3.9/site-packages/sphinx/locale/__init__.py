"""
    sphinx.locale
    ~~~~~~~~~~~~~

    Locale utilities.

    :copyright: Copyright 2007-2022 by the Sphinx team, see AUTHORS.
    :license: BSD, see LICENSE for details.
"""

import gettext
import locale
from collections import UserString, defaultdict
from gettext import NullTranslations
from typing import Any, Callable, Dict, Iterable, List, Optional, Tuple, Union


class _TranslationProxy(UserString):
    """
    Class for proxy strings from gettext translations. This is a helper for the
    lazy_* functions from this module.

    The proxy implementation attempts to be as complete as possible, so that
    the lazy objects should mostly work as expected, for example for sorting.

    This inherits from UserString because some docutils versions use UserString
    for their Text nodes, which then checks its argument for being either a
    basestring or UserString, otherwise calls str() -- not unicode() -- on it.
    """
    __slots__ = ('_func', '_args')

    def __new__(cls, func: Callable, *args: str) -> object:  # type: ignore
        if not args:
            # not called with "function" and "arguments", but a plain string
            return str(func)
        return object.__new__(cls)

    def __getnewargs__(self) -> Tuple[str]:
        return (self._func,) + self._args  # type: ignore

    def __init__(self, func: Callable, *args: str) -> None:
        self._func = func
        self._args = args

    @property
    def data(self) -> str:  # type: ignore
        return self._func(*self._args)

    # replace function from UserString; it instantiates a self.__class__
    # for the encoding result

    def encode(self, encoding: str = None, errors: str = None) -> bytes:  # type: ignore
        if encoding:
            if errors:
                return self.data.encode(encoding, errors)
            else:
                return self.data.encode(encoding)
        else:
            return self.data.encode()

    def __dir__(self) -> List[str]:
        return dir(str)

    def __str__(self) -> str:
        return str(self.data)

    def __add__(self, other: str) -> str:  # type: ignore
        return self.data + other

    def __radd__(self, other: str) -> str:
        return other + self.data

    def __mod__(self, other: str) -> str:  # type: ignore
        return self.data % other

    def __rmod__(self, other: str) -> str:
        return other % self.data

    def __mul__(self, other: Any) -> str:  # type: ignore
        return self.data * other

    def __rmul__(self, other: Any) -> str:
        return other * self.data

    def __getattr__(self, name: str) -> Any:
        if name == '__members__':
            return self.__dir__()
        return getattr(self.data, name)

    def __getstate__(self) -> Tuple[Callable, Tuple[str, ...]]:
        return self._func, self._args

    def __setstate__(self, tup: Tuple[Callable, Tuple[str]]) -> None:
        self._func, self._args = tup

    def __copy__(self) -> "_TranslationProxy":
        return self

    def __repr__(self) -> str:
        try:
            return 'i' + repr(str(self.data))
        except Exception:
            return '<%s broken>' % self.__class__.__name__


translators: Dict[Tuple[str, str], NullTranslations] = defaultdict(NullTranslations)


def init(locale_dirs: List[Optional[str]], language: Optional[str],
         catalog: str = 'sphinx', namespace: str = 'general') -> Tuple[NullTranslations, bool]:
    """Look for message catalogs in `locale_dirs` and *ensure* that there is at
    least a NullTranslations catalog set in `translators`. If called multiple
    times or if several ``.mo`` files are found, their contents are merged
    together (thus making ``init`` reentrant).
    """
    global translators
    translator = translators.get((namespace, catalog))
    # ignore previously failed attempts to find message catalogs
    if translator.__class__ is NullTranslations:
        translator = None
    # the None entry is the system's default locale path
    has_translation = True

    if language and '_' in language:
        # for language having country code (like "de_AT")
        languages: Optional[List[str]] = [language, language.split('_')[0]]
    elif language:
        languages = [language]
    else:
        languages = None

    # loading
    for dir_ in locale_dirs:
        try:
            trans = gettext.translation(catalog, localedir=dir_, languages=languages)
            if translator is None:
                translator = trans
            else:
                translator.add_fallback(trans)
        except Exception:
            # Language couldn't be found in the specified path
            pass
    # guarantee translators[(namespace, catalog)] exists
    if translator is None:
        translator = NullTranslations()
        has_translation = False
    translators[(namespace, catalog)] = translator
    return translator, has_translation


def setlocale(category: int, value: Union[str, Iterable[str]] = None) -> None:
    """Update locale settings.

    This does not throw any exception even if update fails.
    This is workaround for Python's bug.

    For more details:

    * https://github.com/sphinx-doc/sphinx/issues/5724
    * https://bugs.python.org/issue18378#msg215215

    .. note:: Only for internal use.  Please don't call this method from extensions.
              This will be removed in future.
    """
    try:
        locale.setlocale(category, value)
    except locale.Error:
        pass


def init_console(locale_dir: str, catalog: str) -> Tuple[NullTranslations, bool]:
    """Initialize locale for console.

    .. versionadded:: 1.8
    """
    try:
        # encoding is ignored
        language, _ = locale.getlocale(locale.LC_MESSAGES)  # type: Tuple[Optional[str], Any]
    except AttributeError:
        # LC_MESSAGES is not always defined. Fallback to the default language
        # in case it is not.
        language = None
    return init([locale_dir], language, catalog, 'console')


def get_translator(catalog: str = 'sphinx', namespace: str = 'general') -> NullTranslations:
    return translators[(namespace, catalog)]


def is_translator_registered(catalog: str = 'sphinx', namespace: str = 'general') -> bool:
    return (namespace, catalog) in translators


def _lazy_translate(catalog: str, namespace: str, message: str) -> str:
    """Used instead of _ when creating TranslationProxy, because _ is
    not bound yet at that time.
    """
    translator = get_translator(catalog, namespace)
    return translator.gettext(message)


def get_translation(catalog: str, namespace: str = 'general') -> Callable:
    """Get a translation function based on the *catalog* and *namespace*.

    The extension can use this API to translate the messages on the
    extension::

        import os
        from sphinx.locale import get_translation

        MESSAGE_CATALOG_NAME = 'myextension'  # name of *.pot, *.po and *.mo files
        _ = get_translation(MESSAGE_CATALOG_NAME)
        text = _('Hello Sphinx!')


        def setup(app):
            package_dir = os.path.abspath(os.path.dirname(__file__))
            locale_dir = os.path.join(package_dir, 'locales')
            app.add_message_catalog(MESSAGE_CATALOG_NAME, locale_dir)

    With this code, sphinx searches a message catalog from
    ``${package_dir}/locales/${language}/LC_MESSAGES/myextension.mo``.
    The :confval:`language` is used for the searching.

    .. versionadded:: 1.8
    """
    def gettext(message: str, *args: Any) -> str:
        if not is_translator_registered(catalog, namespace):
            # not initialized yet
            return _TranslationProxy(_lazy_translate, catalog, namespace, message)  # type: ignore  # NOQA
        else:
            translator = get_translator(catalog, namespace)
            if len(args) <= 1:
                return translator.gettext(message)
            else:  # support pluralization
                return translator.ngettext(message, args[0], args[1])

    return gettext


# A shortcut for sphinx-core
#: Translation function for messages on documentation (menu, labels, themes and so on).
#: This function follows :confval:`language` setting.
_ = get_translation('sphinx')
#: Translation function for console messages
#: This function follows locale setting (`LC_ALL`, `LC_MESSAGES` and so on).
__ = get_translation('sphinx', 'console')


# labels
admonitionlabels = {
    'attention': _('Attention'),
    'caution':   _('Caution'),
    'danger':    _('Danger'),
    'error':     _('Error'),
    'hint':      _('Hint'),
    'important': _('Important'),
    'note':      _('Note'),
    'seealso':   _('See also'),
    'tip':       _('Tip'),
    'warning':   _('Warning'),
}

# Moved to sphinx.directives.other (will be overridden later)
versionlabels: Dict[str, str] = {}

# Moved to sphinx.domains.python (will be overridden later)
pairindextypes: Dict[str, str] = {}
