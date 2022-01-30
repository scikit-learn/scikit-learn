import argparse
import configparser
import glob as fileglob
from io import StringIO
import os
import re
import sys

import tomli
from typing import (Any, Callable, Dict, List, Mapping, MutableMapping,  Optional, Sequence,
                    TextIO, Tuple, Union)
from typing_extensions import Final, TypeAlias as _TypeAlias

from mypy import defaults
from mypy.options import Options, PER_MODULE_OPTIONS

_CONFIG_VALUE_TYPES: _TypeAlias = Union[
    str, bool, int, float, Dict[str, str], List[str], Tuple[int, int],
]
_INI_PARSER_CALLABLE: _TypeAlias = Callable[[Any], _CONFIG_VALUE_TYPES]


def parse_version(v: str) -> Tuple[int, int]:
    m = re.match(r'\A(\d)\.(\d+)\Z', v)
    if not m:
        raise argparse.ArgumentTypeError(
            "Invalid python version '{}' (expected format: 'x.y')".format(v))
    major, minor = int(m.group(1)), int(m.group(2))
    if major == 2:
        if minor != 7:
            raise argparse.ArgumentTypeError(
                "Python 2.{} is not supported (must be 2.7)".format(minor))
    elif major == 3:
        if minor < defaults.PYTHON3_VERSION_MIN[1]:
            raise argparse.ArgumentTypeError(
                "Python 3.{0} is not supported (must be {1}.{2} or higher)".format(minor,
                                                                    *defaults.PYTHON3_VERSION_MIN))
    else:
        raise argparse.ArgumentTypeError(
            "Python major version '{}' out of range (must be 2 or 3)".format(major))
    return major, minor


def try_split(v: Union[str, Sequence[str]], split_regex: str = '[,]') -> List[str]:
    """Split and trim a str or list of str into a list of str"""
    if isinstance(v, str):
        return [p.strip() for p in re.split(split_regex, v)]

    return [p.strip() for p in v]


def expand_path(path: str) -> str:
    """Expand the user home directory and any environment variables contained within
    the provided path.
    """

    return os.path.expandvars(os.path.expanduser(path))


def str_or_array_as_list(v: Union[str, Sequence[str]]) -> List[str]:
    if isinstance(v, str):
        return [v.strip()] if v.strip() else []
    return [p.strip() for p in v if p.strip()]


def split_and_match_files_list(paths: Sequence[str]) -> List[str]:
    """Take a list of files/directories (with support for globbing through the glob library).

    Where a path/glob matches no file, we still include the raw path in the resulting list.

    Returns a list of file paths
    """
    expanded_paths = []

    for path in paths:
        path = expand_path(path.strip())
        globbed_files = fileglob.glob(path, recursive=True)
        if globbed_files:
            expanded_paths.extend(globbed_files)
        else:
            expanded_paths.append(path)

    return expanded_paths


def split_and_match_files(paths: str) -> List[str]:
    """Take a string representing a list of files/directories (with support for globbing
    through the glob library).

    Where a path/glob matches no file, we still include the raw path in the resulting list.

    Returns a list of file paths
    """

    return split_and_match_files_list(paths.split(','))


def check_follow_imports(choice: str) -> str:
    choices = ['normal', 'silent', 'skip', 'error']
    if choice not in choices:
        raise argparse.ArgumentTypeError(
            "invalid choice '{}' (choose from {})".format(
                choice,
                ', '.join("'{}'".format(x) for x in choices)))
    return choice


# For most options, the type of the default value set in options.py is
# sufficient, and we don't have to do anything here.  This table
# exists to specify types for values initialized to None or container
# types.
ini_config_types: Final[Dict[str, _INI_PARSER_CALLABLE]] = {
    'python_version': parse_version,
    'strict_optional_whitelist': lambda s: s.split(),
    'custom_typing_module': str,
    'custom_typeshed_dir': expand_path,
    'mypy_path': lambda s: [expand_path(p.strip()) for p in re.split('[,:]', s)],
    'files': split_and_match_files,
    'quickstart_file': expand_path,
    'junit_xml': expand_path,
    # These two are for backwards compatibility
    'silent_imports': bool,
    'almost_silent': bool,
    'follow_imports': check_follow_imports,
    'no_site_packages': bool,
    'plugins': lambda s: [p.strip() for p in s.split(',')],
    'always_true': lambda s: [p.strip() for p in s.split(',')],
    'always_false': lambda s: [p.strip() for p in s.split(',')],
    'disable_error_code': lambda s: [p.strip() for p in s.split(',')],
    'enable_error_code': lambda s: [p.strip() for p in s.split(',')],
    'package_root': lambda s: [p.strip() for p in s.split(',')],
    'cache_dir': expand_path,
    'python_executable': expand_path,
    'strict': bool,
    'exclude': lambda s: [s.strip()],
}

# Reuse the ini_config_types and overwrite the diff
toml_config_types: Final[Dict[str, _INI_PARSER_CALLABLE]] = ini_config_types.copy()
toml_config_types.update({
    'python_version': lambda s: parse_version(str(s)),
    'strict_optional_whitelist': try_split,
    'mypy_path': lambda s: [expand_path(p) for p in try_split(s, '[,:]')],
    'files': lambda s: split_and_match_files_list(try_split(s)),
    'follow_imports': lambda s: check_follow_imports(str(s)),
    'plugins': try_split,
    'always_true': try_split,
    'always_false': try_split,
    'disable_error_code': try_split,
    'enable_error_code': try_split,
    'package_root': try_split,
    'exclude': str_or_array_as_list,
})


def parse_config_file(options: Options, set_strict_flags: Callable[[], None],
                      filename: Optional[str],
                      stdout: Optional[TextIO] = None,
                      stderr: Optional[TextIO] = None) -> None:
    """Parse a config file into an Options object.

    Errors are written to stderr but are not fatal.

    If filename is None, fall back to default config files.
    """
    stdout = stdout or sys.stdout
    stderr = stderr or sys.stderr

    if filename is not None:
        config_files: Tuple[str, ...] = (filename,)
    else:
        config_files = tuple(map(os.path.expanduser, defaults.CONFIG_FILES))

    config_parser = configparser.RawConfigParser()

    for config_file in config_files:
        if not os.path.exists(config_file):
            continue
        try:
            if is_toml(config_file):
                with open(config_file, encoding="utf-8") as f:
                    toml_data = tomli.loads(f.read())
                # Filter down to just mypy relevant toml keys
                toml_data = toml_data.get('tool', {})
                if 'mypy' not in toml_data:
                    continue
                toml_data = {'mypy': toml_data['mypy']}
                parser: MutableMapping[str, Any] = destructure_overrides(toml_data)
                config_types = toml_config_types
            else:
                config_parser.read(config_file)
                parser = config_parser
                config_types = ini_config_types
        except (tomli.TOMLDecodeError, configparser.Error, ConfigTOMLValueError) as err:
            print("%s: %s" % (config_file, err), file=stderr)
        else:
            if config_file in defaults.SHARED_CONFIG_FILES and 'mypy' not in parser:
                continue
            file_read = config_file
            options.config_file = file_read
            break
    else:
        return

    os.environ['MYPY_CONFIG_FILE_DIR'] = os.path.dirname(
            os.path.abspath(config_file))

    if 'mypy' not in parser:
        if filename or file_read not in defaults.SHARED_CONFIG_FILES:
            print("%s: No [mypy] section in config file" % file_read, file=stderr)
    else:
        section = parser['mypy']
        prefix = '%s: [%s]: ' % (file_read, 'mypy')
        updates, report_dirs = parse_section(
            prefix, options, set_strict_flags, section, config_types, stderr)
        for k, v in updates.items():
            setattr(options, k, v)
        options.report_dirs.update(report_dirs)

    for name, section in parser.items():
        if name.startswith('mypy-'):
            prefix = get_prefix(file_read, name)
            updates, report_dirs = parse_section(
                prefix, options, set_strict_flags, section, config_types, stderr)
            if report_dirs:
                print("%sPer-module sections should not specify reports (%s)" %
                      (prefix, ', '.join(s + '_report' for s in sorted(report_dirs))),
                      file=stderr)
            if set(updates) - PER_MODULE_OPTIONS:
                print("%sPer-module sections should only specify per-module flags (%s)" %
                      (prefix, ', '.join(sorted(set(updates) - PER_MODULE_OPTIONS))),
                      file=stderr)
                updates = {k: v for k, v in updates.items() if k in PER_MODULE_OPTIONS}
            globs = name[5:]
            for glob in globs.split(','):
                # For backwards compatibility, replace (back)slashes with dots.
                glob = glob.replace(os.sep, '.')
                if os.altsep:
                    glob = glob.replace(os.altsep, '.')

                if (any(c in glob for c in '?[]!') or
                        any('*' in x and x != '*' for x in glob.split('.'))):
                    print("%sPatterns must be fully-qualified module names, optionally "
                          "with '*' in some components (e.g spam.*.eggs.*)"
                          % prefix,
                          file=stderr)
                else:
                    options.per_module_options[glob] = updates


def get_prefix(file_read: str, name: str) -> str:
    if is_toml(file_read):
        module_name_str = 'module = "%s"' % '-'.join(name.split('-')[1:])
    else:
        module_name_str = name

    return '%s: [%s]: ' % (file_read, module_name_str)


def is_toml(filename: str) -> bool:
    return filename.lower().endswith('.toml')


def destructure_overrides(toml_data: Dict[str, Any]) -> Dict[str, Any]:
    """Take the new [[tool.mypy.overrides]] section array in the pyproject.toml file,
    and convert it back to a flatter structure that the existing config_parser can handle.

    E.g. the following pyproject.toml file:

        [[tool.mypy.overrides]]
        module = [
            "a.b",
            "b.*"
        ]
        disallow_untyped_defs = true

        [[tool.mypy.overrides]]
        module = 'c'
        disallow_untyped_defs = false

    Would map to the following config dict that it would have gotten from parsing an equivalent
    ini file:

        {
            "mypy-a.b": {
                disallow_untyped_defs = true,
            },
            "mypy-b.*": {
                disallow_untyped_defs = true,
            },
            "mypy-c": {
                disallow_untyped_defs: false,
            },
        }
    """
    if 'overrides' not in toml_data['mypy']:
        return toml_data

    if not isinstance(toml_data['mypy']['overrides'], list):
        raise ConfigTOMLValueError("tool.mypy.overrides sections must be an array. Please make "
                         "sure you are using double brackets like so: [[tool.mypy.overrides]]")

    result = toml_data.copy()
    for override in result['mypy']['overrides']:
        if 'module' not in override:
            raise ConfigTOMLValueError("toml config file contains a [[tool.mypy.overrides]] "
                             "section, but no module to override was specified.")

        if isinstance(override['module'], str):
            modules = [override['module']]
        elif isinstance(override['module'], list):
            modules = override['module']
        else:
            raise ConfigTOMLValueError("toml config file contains a [[tool.mypy.overrides]] "
                             "section with a module value that is not a string or a list of "
                             "strings")

        for module in modules:
            module_overrides = override.copy()
            del module_overrides['module']
            old_config_name = 'mypy-%s' % module
            if old_config_name not in result:
                result[old_config_name] = module_overrides
            else:
                for new_key, new_value in module_overrides.items():
                    if (new_key in result[old_config_name] and
                            result[old_config_name][new_key] != new_value):
                        raise ConfigTOMLValueError("toml config file contains "
                                         "[[tool.mypy.overrides]] sections with conflicting "
                                         "values. Module '%s' has two different values for '%s'"
                                         % (module, new_key))
                    result[old_config_name][new_key] = new_value

    del result['mypy']['overrides']
    return result


def parse_section(prefix: str, template: Options,
                  set_strict_flags: Callable[[], None],
                  section: Mapping[str, Any],
                  config_types: Dict[str, Any],
                  stderr: TextIO = sys.stderr
                  ) -> Tuple[Dict[str, object], Dict[str, str]]:
    """Parse one section of a config file.

    Returns a dict of option values encountered, and a dict of report directories.
    """
    results: Dict[str, object] = {}
    report_dirs: Dict[str, str] = {}
    for key in section:
        invert = False
        options_key = key
        if key in config_types:
            ct = config_types[key]
        else:
            dv = None
            # We have to keep new_semantic_analyzer in Options
            # for plugin compatibility but it is not a valid option anymore.
            assert hasattr(template, 'new_semantic_analyzer')
            if key != 'new_semantic_analyzer':
                dv = getattr(template, key, None)
            if dv is None:
                if key.endswith('_report'):
                    report_type = key[:-7].replace('_', '-')
                    if report_type in defaults.REPORTER_NAMES:
                        report_dirs[report_type] = str(section[key])
                    else:
                        print("%sUnrecognized report type: %s" % (prefix, key),
                              file=stderr)
                    continue
                if key.startswith('x_'):
                    pass  # Don't complain about `x_blah` flags
                elif key.startswith('no_') and hasattr(template, key[3:]):
                    options_key = key[3:]
                    invert = True
                elif key.startswith('allow') and hasattr(template, 'dis' + key):
                    options_key = 'dis' + key
                    invert = True
                elif key.startswith('disallow') and hasattr(template, key[3:]):
                    options_key = key[3:]
                    invert = True
                elif key == 'strict':
                    pass  # Special handling below
                else:
                    print("%sUnrecognized option: %s = %s" % (prefix, key, section[key]),
                          file=stderr)
                if invert:
                    dv = getattr(template, options_key, None)
                else:
                    continue
            ct = type(dv)
        v: Any = None
        try:
            if ct is bool:
                if isinstance(section, dict):
                    v = convert_to_boolean(section.get(key))
                else:
                    v = section.getboolean(key)  # type: ignore[attr-defined]  # Until better stub
                if invert:
                    v = not v
            elif callable(ct):
                if invert:
                    print("%sCan not invert non-boolean key %s" % (prefix, options_key),
                          file=stderr)
                    continue
                try:
                    v = ct(section.get(key))
                except argparse.ArgumentTypeError as err:
                    print("%s%s: %s" % (prefix, key, err), file=stderr)
                    continue
            else:
                print("%sDon't know what type %s should have" % (prefix, key), file=stderr)
                continue
        except ValueError as err:
            print("%s%s: %s" % (prefix, key, err), file=stderr)
            continue
        if key == 'strict':
            if v:
                set_strict_flags()
            continue
        if key == 'silent_imports':
            print("%ssilent_imports has been replaced by "
                  "ignore_missing_imports=True; follow_imports=skip" % prefix, file=stderr)
            if v:
                if 'ignore_missing_imports' not in results:
                    results['ignore_missing_imports'] = True
                if 'follow_imports' not in results:
                    results['follow_imports'] = 'skip'
        if key == 'almost_silent':
            print("%salmost_silent has been replaced by "
                  "follow_imports=error" % prefix, file=stderr)
            if v:
                if 'follow_imports' not in results:
                    results['follow_imports'] = 'error'
        results[options_key] = v
    return results, report_dirs


def convert_to_boolean(value: Optional[Any]) -> bool:
    """Return a boolean value translating from other types if necessary."""
    if isinstance(value, bool):
        return value
    if not isinstance(value, str):
        value = str(value)
    if value.lower() not in configparser.RawConfigParser.BOOLEAN_STATES:
        raise ValueError('Not a boolean: %s' % value)
    return configparser.RawConfigParser.BOOLEAN_STATES[value.lower()]


def split_directive(s: str) -> Tuple[List[str], List[str]]:
    """Split s on commas, except during quoted sections.

    Returns the parts and a list of error messages."""
    parts = []
    cur: List[str] = []
    errors = []
    i = 0
    while i < len(s):
        if s[i] == ',':
            parts.append(''.join(cur).strip())
            cur = []
        elif s[i] == '"':
            i += 1
            while i < len(s) and s[i] != '"':
                cur.append(s[i])
                i += 1
            if i == len(s):
                errors.append("Unterminated quote in configuration comment")
                cur.clear()
        else:
            cur.append(s[i])
        i += 1
    if cur:
        parts.append(''.join(cur).strip())

    return parts, errors


def mypy_comments_to_config_map(line: str,
                                template: Options) -> Tuple[Dict[str, str], List[str]]:
    """Rewrite the mypy comment syntax into ini file syntax.

    Returns
    """
    options = {}
    entries, errors = split_directive(line)
    for entry in entries:
        if '=' not in entry:
            name = entry
            value = None
        else:
            name, value = [x.strip() for x in entry.split('=', 1)]

        name = name.replace('-', '_')
        if value is None:
            value = 'True'
        options[name] = value

    return options, errors


def parse_mypy_comments(
        args: List[Tuple[int, str]],
        template: Options) -> Tuple[Dict[str, object], List[Tuple[int, str]]]:
    """Parse a collection of inline mypy: configuration comments.

    Returns a dictionary of options to be applied and a list of error messages
    generated.
    """

    errors: List[Tuple[int, str]] = []
    sections = {}

    for lineno, line in args:
        # In order to easily match the behavior for bools, we abuse configparser.
        # Oddly, the only way to get the SectionProxy object with the getboolean
        # method is to create a config parser.
        parser = configparser.RawConfigParser()
        options, parse_errors = mypy_comments_to_config_map(line, template)
        parser['dummy'] = options
        errors.extend((lineno, x) for x in parse_errors)

        stderr = StringIO()
        strict_found = False

        def set_strict_flags() -> None:
            nonlocal strict_found
            strict_found = True

        new_sections, reports = parse_section(
            '', template, set_strict_flags, parser['dummy'], ini_config_types, stderr=stderr)
        errors.extend((lineno, x) for x in stderr.getvalue().strip().split('\n') if x)
        if reports:
            errors.append((lineno, "Reports not supported in inline configuration"))
        if strict_found:
            errors.append((lineno,
                           'Setting "strict" not supported in inline configuration: specify it in '
                           'a configuration file instead, or set individual inline flags '
                           '(see "mypy -h" for the list of flags enabled in strict mode)'))

        sections.update(new_sections)

    return sections, errors


def get_config_module_names(filename: Optional[str], modules: List[str]) -> str:
    if not filename or not modules:
        return ''

    if not is_toml(filename):
        return ", ".join("[mypy-%s]" % module for module in modules)

    return "module = ['%s']" % ("', '".join(sorted(modules)))


class ConfigTOMLValueError(ValueError):
    pass
