"""
    sphinx.theming
    ~~~~~~~~~~~~~~

    Theming support for HTML builders.

    :copyright: Copyright 2007-2022 by the Sphinx team, see AUTHORS.
    :license: BSD, see LICENSE for details.
"""

import configparser
import os
import shutil
import tempfile
from os import path
from typing import TYPE_CHECKING, Any, Dict, List
from zipfile import ZipFile

try:  # Python < 3.10 (backport)
    from importlib_metadata import entry_points
except ImportError:
    from importlib.metadata import entry_points

from sphinx import package_dir
from sphinx.errors import ThemeError
from sphinx.locale import __
from sphinx.util import logging
from sphinx.util.osutil import ensuredir

if TYPE_CHECKING:
    from sphinx.application import Sphinx


logger = logging.getLogger(__name__)

NODEFAULT = object()
THEMECONF = 'theme.conf'


def extract_zip(filename: str, targetdir: str) -> None:
    """Extract zip file to target directory."""
    ensuredir(targetdir)

    with ZipFile(filename) as archive:
        for name in archive.namelist():
            if name.endswith('/'):
                continue
            entry = path.join(targetdir, name)
            ensuredir(path.dirname(entry))
            with open(path.join(entry), 'wb') as fp:
                fp.write(archive.read(name))


class Theme:
    """A Theme is a set of HTML templates and configurations.

    This class supports both theme directory and theme archive (zipped theme)."""

    def __init__(self, name: str, theme_path: str, factory: "HTMLThemeFactory") -> None:
        self.name = name
        self.base = None
        self.rootdir = None

        if path.isdir(theme_path):
            # already a directory, do nothing
            self.rootdir = None
            self.themedir = theme_path
        else:
            # extract the theme to a temp directory
            self.rootdir = tempfile.mkdtemp('sxt')
            self.themedir = path.join(self.rootdir, name)
            extract_zip(theme_path, self.themedir)

        self.config = configparser.RawConfigParser()
        self.config.read(path.join(self.themedir, THEMECONF))

        try:
            inherit = self.config.get('theme', 'inherit')
        except configparser.NoSectionError as exc:
            raise ThemeError(__('theme %r doesn\'t have "theme" setting') % name) from exc
        except configparser.NoOptionError as exc:
            raise ThemeError(__('theme %r doesn\'t have "inherit" setting') % name) from exc

        if inherit != 'none':
            try:
                self.base = factory.create(inherit)
            except ThemeError as exc:
                raise ThemeError(__('no theme named %r found, inherited by %r') %
                                 (inherit, name)) from exc

    def get_theme_dirs(self) -> List[str]:
        """Return a list of theme directories, beginning with this theme's,
        then the base theme's, then that one's base theme's, etc.
        """
        if self.base is None:
            return [self.themedir]
        else:
            return [self.themedir] + self.base.get_theme_dirs()

    def get_config(self, section: str, name: str, default: Any = NODEFAULT) -> Any:
        """Return the value for a theme configuration setting, searching the
        base theme chain.
        """
        try:
            return self.config.get(section, name)
        except (configparser.NoOptionError, configparser.NoSectionError) as exc:
            if self.base:
                return self.base.get_config(section, name, default)

            if default is NODEFAULT:
                raise ThemeError(__('setting %s.%s occurs in none of the '
                                    'searched theme configs') % (section, name)) from exc
            else:
                return default

    def get_options(self, overrides: Dict[str, Any] = {}) -> Dict[str, Any]:
        """Return a dictionary of theme options and their values."""
        if self.base:
            options = self.base.get_options()
        else:
            options = {}

        try:
            options.update(self.config.items('options'))
        except configparser.NoSectionError:
            pass

        for option, value in overrides.items():
            if option not in options:
                logger.warning(__('unsupported theme option %r given') % option)
            else:
                options[option] = value

        return options

    def cleanup(self) -> None:
        """Remove temporary directories."""
        if self.rootdir:
            try:
                shutil.rmtree(self.rootdir)
            except Exception:
                pass
        if self.base:
            self.base.cleanup()


def is_archived_theme(filename: str) -> bool:
    """Check whether the specified file is an archived theme file or not."""
    try:
        with ZipFile(filename) as f:
            return THEMECONF in f.namelist()
    except Exception:
        return False


class HTMLThemeFactory:
    """A factory class for HTML Themes."""

    def __init__(self, app: "Sphinx") -> None:
        self.app = app
        self.themes = app.registry.html_themes
        self.load_builtin_themes()
        if getattr(app.config, 'html_theme_path', None):
            self.load_additional_themes(app.config.html_theme_path)

    def load_builtin_themes(self) -> None:
        """Load built-in themes."""
        themes = self.find_themes(path.join(package_dir, 'themes'))
        for name, theme in themes.items():
            self.themes[name] = theme

    def load_additional_themes(self, theme_paths: str) -> None:
        """Load additional themes placed at specified directories."""
        for theme_path in theme_paths:
            abs_theme_path = path.abspath(path.join(self.app.confdir, theme_path))
            themes = self.find_themes(abs_theme_path)
            for name, theme in themes.items():
                self.themes[name] = theme

    def load_extra_theme(self, name: str) -> None:
        """Try to load a theme with the specified name."""
        if name == 'alabaster':
            self.load_alabaster_theme()
        else:
            self.load_external_theme(name)

    def load_alabaster_theme(self) -> None:
        """Load alabaster theme."""
        import alabaster
        self.themes['alabaster'] = path.join(alabaster.get_path(), 'alabaster')

    def load_sphinx_rtd_theme(self) -> None:
        """Load sphinx_rtd_theme theme (if installed)."""
        try:
            import sphinx_rtd_theme
            theme_path = sphinx_rtd_theme.get_html_theme_path()
            self.themes['sphinx_rtd_theme'] = path.join(theme_path, 'sphinx_rtd_theme')
        except ImportError:
            pass

    def load_external_theme(self, name: str) -> None:
        """Try to load a theme using entry_points.

        Sphinx refers to ``sphinx_themes`` entry_points.
        """
        # look up for new styled entry_points at first
        theme_entry_points = entry_points(group='sphinx.html_themes')
        try:
            entry_point = theme_entry_points[name]
            self.app.registry.load_extension(self.app, entry_point.module)
            self.app.config.post_init_values()
            return
        except KeyError:
            pass

    def find_themes(self, theme_path: str) -> Dict[str, str]:
        """Search themes from specified directory."""
        themes: Dict[str, str] = {}
        if not path.isdir(theme_path):
            return themes

        for entry in os.listdir(theme_path):
            pathname = path.join(theme_path, entry)
            if path.isfile(pathname) and entry.lower().endswith('.zip'):
                if is_archived_theme(pathname):
                    name = entry[:-4]
                    themes[name] = pathname
                else:
                    logger.warning(__('file %r on theme path is not a valid '
                                      'zipfile or contains no theme'), entry)
            else:
                if path.isfile(path.join(pathname, THEMECONF)):
                    themes[entry] = pathname

        return themes

    def create(self, name: str) -> Theme:
        """Create an instance of theme."""
        if name not in self.themes:
            self.load_extra_theme(name)

        if name not in self.themes and name == 'sphinx_rtd_theme':
            # sphinx_rtd_theme (< 0.2.5)  # RemovedInSphinx60Warning
            logger.warning(__('sphinx_rtd_theme (< 0.3.0) found. '
                              'It will not be available since Sphinx-6.0'))
            self.load_sphinx_rtd_theme()

        if name not in self.themes:
            raise ThemeError(__('no theme named %r found (missing theme.conf?)') % name)

        return Theme(name, self.themes[name], factory=self)
