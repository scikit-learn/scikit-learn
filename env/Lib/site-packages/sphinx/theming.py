"""Theming support for HTML builders."""

from __future__ import annotations

__all__ = ('Theme', 'HTMLThemeFactory')

import configparser
import contextlib
import os
import shutil
import sys
import tempfile
from importlib.metadata import entry_points
from os import path
from typing import TYPE_CHECKING, Any
from zipfile import ZipFile

from sphinx import package_dir
from sphinx.config import check_confval_types as _config_post_init
from sphinx.errors import ThemeError
from sphinx.locale import __
from sphinx.util import logging
from sphinx.util.osutil import ensuredir

if sys.version_info >= (3, 11):
    import tomllib
else:
    import tomli as tomllib


if TYPE_CHECKING:
    from collections.abc import Callable
    from typing import TypedDict

    from typing_extensions import Required

    from sphinx.application import Sphinx

    class _ThemeToml(TypedDict, total=False):
        theme: Required[_ThemeTomlTheme]
        options: dict[str, str]

    class _ThemeTomlTheme(TypedDict, total=False):
        inherit: Required[str]
        stylesheets: list[str]
        sidebars: list[str]
        pygments_style: _ThemeTomlThemePygments

    class _ThemeTomlThemePygments(TypedDict, total=False):
        default: str
        dark: str


logger = logging.getLogger(__name__)

_NO_DEFAULT = object()
_THEME_TOML = 'theme.toml'
_THEME_CONF = 'theme.conf'


class Theme:
    """A Theme is a set of HTML templates and configurations.

    This class supports both theme directory and theme archive (zipped theme).
    """

    def __init__(
        self,
        name: str,
        *,
        configs: dict[str, _ConfigFile],
        paths: list[str],
        tmp_dirs: list[str],
    ) -> None:
        self.name = name
        self._dirs = tuple(paths)
        self._tmp_dirs = tmp_dirs

        options: dict[str, Any] = {}
        self.stylesheets: tuple[str, ...] = ()
        self.sidebar_templates: tuple[str, ...] = ()
        self.pygments_style_default: str | None = None
        self.pygments_style_dark: str | None = None
        for config in reversed(configs.values()):
            options |= config.options
            if config.stylesheets is not None:
                self.stylesheets = config.stylesheets
            if config.sidebar_templates is not None:
                self.sidebar_templates = config.sidebar_templates
            if config.pygments_style_default is not None:
                self.pygments_style_default = config.pygments_style_default
            if config.pygments_style_dark is not None:
                self.pygments_style_dark = config.pygments_style_dark

        self._options = options

    def get_theme_dirs(self) -> list[str]:
        """Return a list of theme directories, beginning with this theme's,
        then the base theme's, then that one's base theme's, etc.
        """
        return list(self._dirs)

    def get_config(self, section: str, name: str, default: Any = _NO_DEFAULT) -> Any:
        """Return the value for a theme configuration setting, searching the
        base theme chain.
        """
        if section == 'theme':
            if name == 'stylesheet':
                value = ', '.join(self.stylesheets) or default
            elif name == 'sidebars':
                value = ', '.join(self.sidebar_templates) or default
            elif name == 'pygments_style':
                value = self.pygments_style_default or default
            elif name == 'pygments_dark_style':
                value = self.pygments_style_dark or default
            else:
                value = default
        elif section == 'options':
            value = self._options.get(name, default)
        else:
            msg = __(
                'Theme configuration sections other than [theme] and [options] '
                'are not supported (tried to get a value from %r).'
            )
            raise ThemeError(msg)
        if value is _NO_DEFAULT:
            msg = __('setting %s.%s occurs in none of the searched theme configs') % (
                section,
                name,
            )
            raise ThemeError(msg)
        return value

    def get_options(self, overrides: dict[str, Any] | None = None) -> dict[str, Any]:
        """Return a dictionary of theme options and their values."""
        if overrides is None:
            overrides = {}

        options = self._options.copy()
        for option, value in overrides.items():
            if option not in options:
                logger.warning(__('unsupported theme option %r given'), option)
            else:
                options[option] = value

        return options

    def _cleanup(self) -> None:
        """Remove temporary directories."""
        for tmp_dir in self._tmp_dirs:
            with contextlib.suppress(Exception):
                shutil.rmtree(tmp_dir)


class HTMLThemeFactory:
    """A factory class for HTML Themes."""

    def __init__(self, app: Sphinx) -> None:
        self._app = app
        self._themes = app.registry.html_themes
        self._entry_point_themes: dict[str, Callable[[], None]] = {}
        self._load_builtin_themes()
        if getattr(app.config, 'html_theme_path', None):
            self._load_additional_themes(app.config.html_theme_path)
        self._load_entry_point_themes()

    def _load_builtin_themes(self) -> None:
        """Load built-in themes."""
        themes = self._find_themes(path.join(package_dir, 'themes'))
        for name, theme in themes.items():
            self._themes[name] = theme

    def _load_additional_themes(self, theme_paths: list[str]) -> None:
        """Load additional themes placed at specified directories."""
        for theme_path in theme_paths:
            abs_theme_path = path.abspath(path.join(self._app.confdir, theme_path))
            themes = self._find_themes(abs_theme_path)
            for name, theme in themes.items():
                self._themes[name] = theme

    def _load_entry_point_themes(self) -> None:
        """Try to load a theme with the specified name.

        This uses the ``sphinx.html_themes`` entry point from package metadata.
        """
        for entry_point in entry_points(group='sphinx.html_themes'):
            if entry_point.name in self._themes:
                continue  # don't overwrite loaded themes

            def _load_theme_closure(
                # bind variables in the function definition
                app: Sphinx = self._app,
                theme_module: str = entry_point.module,
            ) -> None:
                app.setup_extension(theme_module)
                _config_post_init(app, app.config)

            self._entry_point_themes[entry_point.name] = _load_theme_closure

    @staticmethod
    def _find_themes(theme_path: str) -> dict[str, str]:
        """Search themes from specified directory."""
        themes: dict[str, str] = {}
        if not path.isdir(theme_path):
            return themes

        for entry in os.listdir(theme_path):
            pathname = path.join(theme_path, entry)
            if path.isfile(pathname) and entry.lower().endswith('.zip'):
                if _is_archived_theme(pathname):
                    name = entry[:-4]
                    themes[name] = pathname
                else:
                    logger.warning(
                        __(
                            'file %r on theme path is not a valid '
                            'zipfile or contains no theme'
                        ),
                        entry,
                    )
            else:
                toml_path = path.join(pathname, _THEME_TOML)
                conf_path = path.join(pathname, _THEME_CONF)
                if path.isfile(toml_path) or path.isfile(conf_path):
                    themes[entry] = pathname

        return themes

    def create(self, name: str) -> Theme:
        """Create an instance of theme."""
        if name in self._entry_point_themes:
            # Load a deferred theme from an entry point
            entry_point_loader = self._entry_point_themes[name]
            entry_point_loader()
        if name not in self._themes:
            raise ThemeError(__('no theme named %r found (missing theme.toml?)') % name)

        themes, theme_dirs, tmp_dirs = _load_theme_with_ancestors(
            name,
            self._themes,
            self._entry_point_themes,
        )
        return Theme(name, configs=themes, paths=theme_dirs, tmp_dirs=tmp_dirs)


def _is_archived_theme(filename: str, /) -> bool:
    """Check whether the specified file is an archived theme file or not."""
    try:
        with ZipFile(filename) as f:
            namelist = frozenset(f.namelist())
            return _THEME_TOML in namelist or _THEME_CONF in namelist
    except Exception:
        return False


def _load_theme_with_ancestors(
    name: str,
    theme_paths: dict[str, str],
    entry_point_themes: dict[str, Callable[[], None]],
    /,
) -> tuple[dict[str, _ConfigFile], list[str], list[str]]:
    themes: dict[str, _ConfigFile] = {}
    theme_dirs: list[str] = []
    tmp_dirs: list[str] = []

    # having 10+ theme ancestors is ludicrous
    for _ in range(10):
        inherit, theme_dir, tmp_dir, config = _load_theme(name, theme_paths[name])
        theme_dirs.append(theme_dir)
        if tmp_dir is not None:
            tmp_dirs.append(tmp_dir)
        themes[name] = config
        if inherit == 'none':
            break
        if inherit in themes:
            msg = __('The %r theme has circular inheritance') % name
            raise ThemeError(msg)
        if inherit in entry_point_themes and inherit not in theme_paths:
            # Load a deferred theme from an entry point
            entry_point_loader = entry_point_themes[inherit]
            entry_point_loader()
        if inherit not in theme_paths:
            msg = __(
                'The %r theme inherits from %r, which is not a loaded theme. '
                'Loaded themes are: %s'
            ) % (name, inherit, ', '.join(sorted(theme_paths)))
            raise ThemeError(msg)
        name = inherit
    else:
        msg = __('The %r theme has too many ancestors') % name
        raise ThemeError(msg)

    return themes, theme_dirs, tmp_dirs


def _load_theme(
    name: str, theme_path: str, /
) -> tuple[str, str, str | None, _ConfigFile]:
    if path.isdir(theme_path):
        # already a directory, do nothing
        tmp_dir = None
        theme_dir = theme_path
    else:
        # extract the theme to a temp directory
        tmp_dir = tempfile.mkdtemp('sxt')
        theme_dir = path.join(tmp_dir, name)
        _extract_zip(theme_path, theme_dir)

    if path.isfile(toml_path := path.join(theme_dir, _THEME_TOML)):
        _cfg_table = _load_theme_toml(toml_path)
        inherit = _validate_theme_toml(_cfg_table, name)
        config = _convert_theme_toml(_cfg_table)
    elif path.isfile(conf_path := path.join(theme_dir, _THEME_CONF)):
        _cfg_parser = _load_theme_conf(conf_path)
        inherit = _validate_theme_conf(_cfg_parser, name)
        config = _convert_theme_conf(_cfg_parser)
    else:
        raise ThemeError(__('no theme configuration file found in %r') % theme_dir)

    return inherit, theme_dir, tmp_dir, config


def _extract_zip(filename: str, target_dir: str, /) -> None:
    """Extract zip file to target directory."""
    ensuredir(target_dir)

    with ZipFile(filename) as archive:
        for name in archive.namelist():
            if name.endswith('/'):
                continue
            entry = path.join(target_dir, name)
            ensuredir(path.dirname(entry))
            with open(path.join(entry), 'wb') as fp:
                fp.write(archive.read(name))


def _load_theme_toml(config_file_path: str, /) -> _ThemeToml:
    with open(config_file_path, encoding='utf-8') as f:
        config_text = f.read()
    c = tomllib.loads(config_text)
    return {s: c[s] for s in ('theme', 'options') if s in c}  # type: ignore[return-value]


def _validate_theme_toml(cfg: _ThemeToml, name: str) -> str:
    if 'theme' not in cfg:
        msg = __('theme %r doesn\'t have the "theme" table') % name
        raise ThemeError(msg)
    theme = cfg['theme']
    if not isinstance(theme, dict):
        msg = __('The %r theme "[theme]" table is not a table') % name
        raise ThemeError(msg)
    inherit = theme.get('inherit', '')
    if not inherit:
        msg = __('The %r theme must define the "theme.inherit" setting') % name
        raise ThemeError(msg)
    if 'options' in cfg:
        if not isinstance(cfg['options'], dict):
            msg = __('The %r theme "[options]" table is not a table') % name
            raise ThemeError(msg)
    return inherit


def _convert_theme_toml(cfg: _ThemeToml, /) -> _ConfigFile:
    theme = cfg['theme']
    if 'stylesheets' in theme:
        stylesheets: tuple[str, ...] | None = tuple(theme['stylesheets'])
    else:
        stylesheets = None
    if 'sidebars' in theme:
        sidebar_templates: tuple[str, ...] | None = tuple(theme['sidebars'])
    else:
        sidebar_templates = None
    pygments_table = theme.get('pygments_style', {})
    if isinstance(pygments_table, str):
        hint = f'pygments_style = {{ default = "{pygments_table}" }}'
        msg = (
            __('The "theme.pygments_style" setting must be a table. Hint: "%s"') % hint
        )
        raise ThemeError(msg)
    pygments_style_default: str | None = pygments_table.get('default')
    pygments_style_dark: str | None = pygments_table.get('dark')
    return _ConfigFile(
        stylesheets=stylesheets,
        sidebar_templates=sidebar_templates,
        pygments_style_default=pygments_style_default,
        pygments_style_dark=pygments_style_dark,
        options=cfg.get('options', {}),
    )


def _load_theme_conf(config_file_path: str, /) -> configparser.RawConfigParser:
    c = configparser.RawConfigParser()
    c.read(config_file_path, encoding='utf-8')
    return c


def _validate_theme_conf(cfg: configparser.RawConfigParser, name: str) -> str:
    if not cfg.has_section('theme'):
        raise ThemeError(__('theme %r doesn\'t have the "theme" table') % name)
    if inherit := cfg.get('theme', 'inherit', fallback=None):
        return inherit
    msg = __('The %r theme must define the "theme.inherit" setting') % name
    raise ThemeError(msg)


def _convert_theme_conf(cfg: configparser.RawConfigParser, /) -> _ConfigFile:
    if stylesheet := cfg.get('theme', 'stylesheet', fallback=''):
        stylesheets: tuple[str, ...] | None = tuple(
            map(str.strip, stylesheet.split(','))
        )
    else:
        stylesheets = None
    if sidebar := cfg.get('theme', 'sidebars', fallback=''):
        sidebar_templates: tuple[str, ...] | None = tuple(
            map(str.strip, sidebar.split(','))
        )
    else:
        sidebar_templates = None
    pygments_style_default: str | None = cfg.get(
        'theme', 'pygments_style', fallback=None
    )
    pygments_style_dark: str | None = cfg.get(
        'theme', 'pygments_dark_style', fallback=None
    )
    options = dict(cfg.items('options')) if cfg.has_section('options') else {}
    return _ConfigFile(
        stylesheets=stylesheets,
        sidebar_templates=sidebar_templates,
        pygments_style_default=pygments_style_default,
        pygments_style_dark=pygments_style_dark,
        options=options,
    )


class _ConfigFile:
    __slots__ = (
        'stylesheets',
        'sidebar_templates',
        'pygments_style_default',
        'pygments_style_dark',
        'options',
    )

    def __init__(
        self,
        stylesheets: tuple[str, ...] | None,
        sidebar_templates: tuple[str, ...] | None,
        pygments_style_default: str | None,
        pygments_style_dark: str | None,
        options: dict[str, str],
    ) -> None:
        self.stylesheets: tuple[str, ...] | None = stylesheets
        self.sidebar_templates: tuple[str, ...] | None = sidebar_templates
        self.pygments_style_default: str | None = pygments_style_default
        self.pygments_style_dark: str | None = pygments_style_dark
        self.options: dict[str, str] = options.copy()

    def __repr__(self) -> str:
        return (
            f'{self.__class__.__qualname__}('
            f'stylesheets={self.stylesheets!r}, '
            f'sidebar_templates={self.sidebar_templates!r}, '
            f'pygments_style_default={self.pygments_style_default!r}, '
            f'pygments_style_dark={self.pygments_style_dark!r}, '
            f'options={self.options!r})'
        )

    def __eq__(self, other: object) -> bool:
        if isinstance(other, _ConfigFile):
            return (
                self.stylesheets == other.stylesheets
                and self.sidebar_templates == other.sidebar_templates
                and self.pygments_style_default == other.pygments_style_default
                and self.pygments_style_dark == other.pygments_style_dark
                and self.options == other.options
            )
        return NotImplemented

    def __hash__(self) -> int:
        return hash((
            self.__class__.__qualname__,
            self.stylesheets,
            self.sidebar_templates,
            self.pygments_style_default,
            self.pygments_style_dark,
            self.options,
        ))


def _migrate_conf_to_toml(argv: list[str]) -> int:
    if argv[:1] != ['conf_to_toml']:
        raise SystemExit(0)
    argv = argv[1:]
    if len(argv) != 1:
        print('Usage: python -m sphinx.theming conf_to_toml <theme path>')  # NoQA: T201
        raise SystemExit(1)
    theme_dir = path.realpath(argv[0])
    conf_path = path.join(theme_dir, _THEME_CONF)
    if not path.isdir(theme_dir) or not path.isfile(conf_path):
        print(  # NoQA: T201
            f'{theme_dir!r} must be a path to a theme directory containing a "theme.conf" file'
        )
        return 1
    _cfg_parser = _load_theme_conf(conf_path)
    if not _cfg_parser.has_section('theme'):
        print('The "theme" table is missing.')  # NoQA: T201
        return 1
    inherit = _cfg_parser.get('theme', 'inherit', fallback=None)
    if not inherit:
        print('The "theme.inherit" setting is missing.')  # NoQA: T201
        return 1

    toml_lines = [
        '[theme]',
        f'inherit = "{inherit}"',
    ]

    stylesheet = _cfg_parser.get('theme', 'stylesheet', fallback=...)
    if stylesheet == '':
        toml_lines.append('stylesheets = []')
    elif stylesheet is not ...:
        toml_lines.append('stylesheets = [')
        toml_lines.extend(f'    "{s}",' for s in map(str.strip, stylesheet.split(',')))
        toml_lines.append(']')

    sidebar = _cfg_parser.get('theme', 'sidebars', fallback=...)
    if sidebar == '':
        toml_lines.append('sidebars = []')
    elif sidebar is not ...:
        toml_lines.append('sidebars = [')
        toml_lines += [f'    "{s}",' for s in map(str.strip, sidebar.split(','))]
        toml_lines.append(']')

    styles = []
    default = _cfg_parser.get('theme', 'pygments_style', fallback=...)
    if default is not ...:
        styles.append(f'default = "{default}"')
    dark = _cfg_parser.get('theme', 'pygments_dark_style', fallback=...)
    if dark is not ...:
        styles.append(f'dark = "{dark}"')
    if styles:
        toml_lines.append('pygments_style = { ' + ', '.join(styles) + ' }')

    if _cfg_parser.has_section('options'):
        toml_lines.append('')
        toml_lines.append('[options]')
        toml_lines += [
            f'{key} = "{d}"'
            for key, default in _cfg_parser.items('options')
            if (d := default.replace('"', r'\"')) or True
        ]

    toml_path = path.join(theme_dir, _THEME_TOML)
    with open(toml_path, 'w', encoding='utf-8') as f:
        f.write('\n'.join(toml_lines) + '\n')
    print(f'Written converted settings to {toml_path!r}')  # NoQA: T201
    return 0


if __name__ == '__main__':
    raise SystemExit(_migrate_conf_to_toml(sys.argv[1:]))
