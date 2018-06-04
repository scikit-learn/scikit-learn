"""Configuration file parsing and utilities."""

import copy
import itertools
import os
from collections import Set, namedtuple
from re import compile as re


from configparser import RawConfigParser


from .utils import __version__, log
from .violations import ErrorRegistry, conventions


def check_initialized(method):
    """Check that the configuration object was initialized."""
    def _decorator(self, *args, **kwargs):
        if self._arguments is None or self._options is None:
            raise RuntimeError('using an uninitialized configuration')
        return method(self, *args, **kwargs)
    return _decorator


class ConfigurationParser(object):
    """Responsible for parsing configuration from files and CLI.

    There are 2 types of configurations: Run configurations and Check
    configurations.

    Run Configurations:
    ------------------
    Responsible for deciding things that are related to the user interface and
    configuration discovery, e.g. verbosity, debug options, etc.
    All run configurations default to `False` or `None` and are decided only 
    by CLI.

    Check Configurations:
    --------------------
    Configurations that are related to which files and errors will be checked.
    These are configurable in 2 ways: using the CLI, and using configuration
    files.

    Configuration files are nested within the file system, meaning that the
    closer a configuration file is to a checked file, the more relevant it will
    be. For instance, imagine this directory structure:

    A
    +-- tox.ini: sets `select=D100`
    +-- B
        +-- foo.py
        +-- tox.ini: sets `add-ignore=D100`

    Then `foo.py` will not be checked for `D100`.
    The configuration build algorithm is described in `self._get_config`.

    Note: If any of `BASE_ERROR_SELECTION_OPTIONS` was selected in the CLI, all
    configuration files will be ignored and each file will be checked for
    the error codes supplied in the CLI.

    """

    CONFIG_FILE_OPTIONS = ('convention', 'select', 'ignore', 'add-select',
                           'add-ignore', 'match', 'match-dir',
                           'ignore-decorators')
    BASE_ERROR_SELECTION_OPTIONS = ('ignore', 'select', 'convention')

    DEFAULT_MATCH_RE = '(?!test_).*\.py'
    DEFAULT_MATCH_DIR_RE = '[^\.].*'
    DEFAULT_IGNORE_DECORATORS_RE = ''
    DEFAULT_CONVENTION = conventions.pep257

    PROJECT_CONFIG_FILES = (
        'setup.cfg',
        'tox.ini',
        '.pydocstyle',
        '.pydocstyle.ini',
        '.pydocstylerc',
        '.pydocstylerc.ini',
        # The following is deprecated, but remains for backwards compatibility.
        '.pep257',
    )

    POSSIBLE_SECTION_NAMES = ('pydocstyle', 'pep257')

    def __init__(self):
        """Create a configuration parser."""
        self._cache = {}
        self._override_by_cli = None
        self._options = self._arguments = self._run_conf = None
        self._parser = self._create_option_parser()

    # ---------------------------- Public Methods -----------------------------

    def get_default_run_configuration(self):
        """Return a `RunConfiguration` object set with default values."""
        options, _ = self._parse_args([])
        return self._create_run_config(options)

    def parse(self):
        """Parse the configuration.

        If one of `BASE_ERROR_SELECTION_OPTIONS` was selected, overrides all
        error codes to check and disregards any error code related
        configurations from the configuration files.

        """
        self._options, self._arguments = self._parse_args()
        self._arguments = self._arguments or ['.']

        if not self._validate_options(self._options):
            raise IllegalConfiguration()

        self._run_conf = self._create_run_config(self._options)

        config = self._create_check_config(self._options, use_defaults=False)
        self._override_by_cli = config

    @check_initialized
    def get_user_run_configuration(self):
        """Return the run configuration for the script."""
        return self._run_conf

    @check_initialized
    def get_files_to_check(self):
        """Generate files and error codes to check on each one.

        Walk dir trees under `self._arguments` and yield file names
        that `match` under each directory that `match_dir`.
        The method locates the configuration for each file name and yields a
        tuple of (filename, [error_codes]).

        With every discovery of a new configuration file `IllegalConfiguration`
        might be raised.

        """
        def _get_matches(config):
            """Return the `match` and `match_dir` functions for `config`."""
            match_func = re(config.match + '$').match
            match_dir_func = re(config.match_dir + '$').match
            return match_func, match_dir_func

        def _get_ignore_decorators(config):
            """Return the `ignore_decorators` as None or regex."""
            if config.ignore_decorators:  # not None and not ''
                ignore_decorators = re(config.ignore_decorators)
            else:
                ignore_decorators = None
            return ignore_decorators

        for name in self._arguments:
            if os.path.isdir(name):
                for root, dirs, filenames in os.walk(name):
                    config = self._get_config(os.path.abspath(root))
                    match, match_dir = _get_matches(config)
                    ignore_decorators = _get_ignore_decorators(config)

                    # Skip any dirs that do not match match_dir
                    dirs[:] = [dir for dir in dirs if match_dir(dir)]

                    for filename in filenames:
                        if match(filename):
                            full_path = os.path.join(root, filename)
                            yield (full_path, list(config.checked_codes),
                                   ignore_decorators)
            else:
                config = self._get_config(os.path.abspath(name))
                match, _ = _get_matches(config)
                ignore_decorators = _get_ignore_decorators(config)
                if match(name):
                    yield (name, list(config.checked_codes), ignore_decorators)

    # --------------------------- Private Methods -----------------------------

    def _get_config_by_discovery(self, node):
        """Get a configuration for checking `node` by config discovery.
        
        Config discovery happens when no explicit config file is specified. The
        file system is searched for config files starting from the directory
        containing the file being checked, and up until the root directory of
        the project.
        
        See `_get_config` for further details.
        
        """
        path = self._get_node_dir(node)

        if path in self._cache:
            return self._cache[path]

        config_file = self._get_config_file_in_folder(path)

        if config_file is None:
            parent_dir, tail = os.path.split(path)
            if tail:
                # No configuration file, simply take the parent's.
                config = self._get_config(parent_dir)
            else:
                # There's no configuration file and no parent directory.
                # Use the default configuration or the one given in the CLI.
                config = self._create_check_config(self._options)
        else:
            # There's a config file! Read it and merge if necessary.
            options, inherit = self._read_configuration_file(config_file)

            parent_dir, tail = os.path.split(path)
            if tail and inherit:
                # There is a parent dir and we should try to merge.
                parent_config = self._get_config(parent_dir)
                config = self._merge_configuration(parent_config, options)
            else:
                # No need to merge or parent dir does not exist.
                config = self._create_check_config(options)

        return config

    def _get_config(self, node):
        """Get and cache the run configuration for `node`.

        If no configuration exists (not local and not for the parent node),
        returns and caches a default configuration.

        The algorithm:
        -------------
        * If the current directory's configuration exists in
           `self._cache` - return it.
        * If a configuration file does not exist in this directory:
          * If the directory is not a root directory:
            * Cache its configuration as this directory's and return it.
          * Else:
            * Cache a default configuration and return it.
        * Else:
          * Read the configuration file.
          * If a parent directory exists AND the configuration file
            allows inheritance:
            * Read the parent configuration by calling this function with the
              parent directory as `node`.
            * Merge the parent configuration with the current one and
              cache it.
        * If the user has specified one of `BASE_ERROR_SELECTION_OPTIONS` in
          the CLI - return the CLI configuration with the configuration match
          clauses
        * Set the `--add-select` and `--add-ignore` CLI configurations.

        """
        if self._run_conf.config is None:
            log.debug('No config file specified, discovering.')
            config = self._get_config_by_discovery(node)
        else:
            log.debug('Using config file %r', self._run_conf.config)
            if not os.path.exists(self._run_conf.config):
                raise IllegalConfiguration('Configuration file {!r} specified '
                                           'via --config was not found.'
                                           .format(self._run_conf.config))

            if None in self._cache:
                return self._cache[None]
            options, _ = self._read_configuration_file(self._run_conf.config)

            if options is None:
                log.warning('Configuration file does not contain a '
                            'pydocstyle section. Using default configuration.')
                config = self._create_check_config(self._options)
            else:
                config = self._create_check_config(options)

        # Make the CLI always win
        final_config = {}
        for attr in CheckConfiguration._fields:
            cli_val = getattr(self._override_by_cli, attr)
            conf_val = getattr(config, attr)
            final_config[attr] = cli_val if cli_val is not None else conf_val

        config = CheckConfiguration(**final_config)

        self._set_add_options(config.checked_codes, self._options)

        # Handle caching
        if self._run_conf.config is not None:
            self._cache[None] = config
        else:
            self._cache[self._get_node_dir(node)] = config
        return config

    @staticmethod
    def _get_node_dir(node):
        """Return the absolute path of the directory of a filesystem node."""
        path = os.path.abspath(node)
        return path if os.path.isdir(path) else os.path.dirname(path)

    def _read_configuration_file(self, path):
        """Try to read and parse `path` as a configuration file.

        If the configurations were illegal (checked with
        `self._validate_options`), raises `IllegalConfiguration`.

        Returns (options, should_inherit).

        """
        parser = RawConfigParser(inline_comment_prefixes=('#', ';'))
        options = None
        should_inherit = True

        if parser.read(path) and self._get_section_name(parser):
            all_options = self._parser.option_list[:]
            for group in self._parser.option_groups:
                all_options.extend(group.option_list)

            option_list = dict([(o.dest, o.type or o.action)
                                for o in all_options])

            # First, read the default values
            new_options, _ = self._parse_args([])

            # Second, parse the configuration
            section_name = self._get_section_name(parser)
            for opt in parser.options(section_name):
                if opt == 'inherit':
                    should_inherit = parser.getboolean(section_name, opt)
                    continue

                if opt.replace('_', '-') not in self.CONFIG_FILE_OPTIONS:
                    log.warning("Unknown option '{}' ignored".format(opt))
                    continue

                normalized_opt = opt.replace('-', '_')
                opt_type = option_list[normalized_opt]
                if opt_type in ('int', 'count'):
                    value = parser.getint(section_name, opt)
                elif opt_type == 'string':
                    value = parser.get(section_name, opt)
                else:
                    assert opt_type in ('store_true', 'store_false')
                    value = parser.getboolean(section_name, opt)
                setattr(new_options, normalized_opt, value)

            # Third, fix the set-options
            options = self._fix_set_options(new_options)

        if options is not None:
            if not self._validate_options(options):
                raise IllegalConfiguration('in file: {}'.format(path))

        return options, should_inherit

    def _merge_configuration(self, parent_config, child_options):
        """Merge parent config into the child options.

        The migration process requires an `options` object for the child in
        order to distinguish between mutually exclusive codes, add-select and
        add-ignore error codes.

        """
        # Copy the parent error codes so we won't override them
        error_codes = copy.deepcopy(parent_config.checked_codes)
        if self._has_exclusive_option(child_options):
            error_codes = self._get_exclusive_error_codes(child_options)

        self._set_add_options(error_codes, child_options)

        kwargs = dict(checked_codes=error_codes)
        for key in ('match', 'match_dir', 'ignore_decorators'):
            kwargs[key] = \
                getattr(child_options, key) or getattr(parent_config, key)
        return CheckConfiguration(**kwargs)

    def _parse_args(self, args=None, values=None):
        """Parse the options using `self._parser` and reformat the options."""
        options, arguments = self._parser.parse_args(args, values)
        return self._fix_set_options(options), arguments

    @staticmethod
    def _create_run_config(options):
        """Create a `RunConfiguration` object from `options`."""
        values = dict([(opt, getattr(options, opt)) for opt in
                       RunConfiguration._fields])
        return RunConfiguration(**values)

    @classmethod
    def _create_check_config(cls, options, use_defaults=True):
        """Create a `CheckConfiguration` object from `options`.

        If `use_defaults`, any of the match options that are `None` will
        be replaced with their default value and the default convention will be
        set for the checked codes.

        """
        checked_codes = None

        if cls._has_exclusive_option(options) or use_defaults:
            checked_codes = cls._get_checked_errors(options)

        kwargs = dict(checked_codes=checked_codes)
        for key in ('match', 'match_dir', 'ignore_decorators'):
            kwargs[key] = getattr(cls, 'DEFAULT_{0}_RE'.format(key.upper())) \
                if getattr(options, key) is None and use_defaults \
                else getattr(options, key)
        return CheckConfiguration(**kwargs)

    @classmethod
    def _get_section_name(cls, parser):
        """Parse options from relevant section."""
        for section_name in cls.POSSIBLE_SECTION_NAMES:
            if parser.has_section(section_name):
                return section_name

        return None

    @classmethod
    def _get_config_file_in_folder(cls, path):
        """Look for a configuration file in `path`.

        If exists return its full path, otherwise None.

        """
        if os.path.isfile(path):
            path = os.path.dirname(path)

        for fn in cls.PROJECT_CONFIG_FILES:
            config = RawConfigParser()
            full_path = os.path.join(path, fn)
            if config.read(full_path) and cls._get_section_name(config):
                return full_path

    @classmethod
    def _get_exclusive_error_codes(cls, options):
        """Extract the error codes from the selected exclusive option."""
        codes = set(ErrorRegistry.get_error_codes())
        checked_codes = None

        if options.ignore is not None:
            ignored = cls._expand_error_codes(options.ignore)
            checked_codes = codes - ignored
        elif options.select is not None:
            checked_codes = cls._expand_error_codes(options.select)
        elif options.convention is not None:
            checked_codes = getattr(conventions, options.convention)

        # To not override the conventions nor the options - copy them.
        return copy.deepcopy(checked_codes)

    @classmethod
    def _set_add_options(cls, checked_codes, options):
        """Set `checked_codes` by the `add_ignore` or `add_select` options."""
        checked_codes |= cls._expand_error_codes(options.add_select)
        checked_codes -= cls._expand_error_codes(options.add_ignore)

    @staticmethod
    def _expand_error_codes(code_parts):
        """Return an expanded set of error codes to ignore."""
        codes = set(ErrorRegistry.get_error_codes())
        expanded_codes = set()

        try:
            for part in code_parts:
                # Dealing with split-lined configurations; The part might begin
                # with a whitespace due to the newline character.
                part = part.strip()
                if not part:
                    continue

                codes_to_add = {code for code in codes
                                if code.startswith(part)}
                if not codes_to_add:
                    log.warn('Error code passed is not a prefix of any known '
                             'errors: %s', part)
                expanded_codes.update(codes_to_add)
        except TypeError as e:
            raise IllegalConfiguration(e)

        return expanded_codes

    @classmethod
    def _get_checked_errors(cls, options):
        """Extract the codes needed to be checked from `options`."""
        checked_codes = cls._get_exclusive_error_codes(options)
        if checked_codes is None:
            checked_codes = cls.DEFAULT_CONVENTION

        cls._set_add_options(checked_codes, options)

        return checked_codes

    @classmethod
    def _validate_options(cls, options):
        """Validate the mutually exclusive options.

        Return `True` iff only zero or one of `BASE_ERROR_SELECTION_OPTIONS`
        was selected.

        """
        for opt1, opt2 in \
                itertools.permutations(cls.BASE_ERROR_SELECTION_OPTIONS, 2):
            if getattr(options, opt1) and getattr(options, opt2):
                log.error('Cannot pass both {} and {}. They are '
                          'mutually exclusive.'.format(opt1, opt2))
                return False

        if options.convention and options.convention not in conventions:
            log.error("Illegal convention '{}'. Possible conventions: {}"
                      .format(options.convention,
                              ', '.join(conventions.keys())))
            return False
        return True

    @classmethod
    def _has_exclusive_option(cls, options):
        """Return `True` iff one or more exclusive options were selected."""
        return any([getattr(options, opt) is not None for opt in
                    cls.BASE_ERROR_SELECTION_OPTIONS])

    @classmethod
    def _fix_set_options(cls, options):
        """Alter the set options from None/strings to sets in place."""
        optional_set_options = ('ignore', 'select')
        mandatory_set_options = ('add_ignore', 'add_select')

        def _get_set(value_str):
            """Split `value_str` by the delimiter `,` and return a set.

            Removes any occurrences of '' in the set.
            Also expand error code prefixes, to avoid doing this for every
            file.

            """
            return cls._expand_error_codes(set(value_str.split(',')) - {''})

        for opt in optional_set_options:
            value = getattr(options, opt)
            if value is not None:
                setattr(options, opt, _get_set(value))

        for opt in mandatory_set_options:
            value = getattr(options, opt)
            if value is None:
                value = ''

            if not isinstance(value, Set):
                value = _get_set(value)

            setattr(options, opt, value)

        return options

    @classmethod
    def _create_option_parser(cls):
        """Return an option parser to parse the command line arguments."""
        from optparse import OptionParser, OptionGroup

        parser = OptionParser(
            version=__version__,
            usage='Usage: pydocstyle [options] [<file|dir>...]')

        option = parser.add_option

        # Run configuration options
        option('-e', '--explain', action='store_true', default=False,
               help='show explanation of each error')
        option('-s', '--source', action='store_true', default=False,
               help='show source for each error')
        option('-d', '--debug', action='store_true', default=False,
               help='print debug information')
        option('-v', '--verbose', action='store_true', default=False,
               help='print status information')
        option('--count', action='store_true', default=False,
               help='print total number of errors to stdout')
        option('--config', metavar='<path>', default=None,
               help='use given config file and disable config discovery')

        check_group = OptionGroup(
            parser,
            'Error Check Options',
            'Only one of --select, --ignore or --convention can be '
            'specified. If none is specified, defaults to '
            '`--convention=pep257`. These three options select the "basic '
            'list" of error codes to check. If you wish to change that list '
            '(for example, if you selected a known convention but wish to '
            'ignore a specific error from it or add a new one) you can '
            'use `--add-[ignore/select]` in order to do so.')
        add_check = check_group.add_option

        # Error check options
        add_check('--select', metavar='<codes>', default=None,
                  help='choose the basic list of checked errors by '
                       'specifying which errors to check for (with a list of '
                       'comma-separated error codes or prefixes). '
                       'for example: --select=D101,D2')
        add_check('--ignore', metavar='<codes>', default=None,
                  help='choose the basic list of checked errors by '
                       'specifying which errors to ignore out of all of the '
                       'available error codes (with a list of '
                       'comma-separated error codes or prefixes). '
                       'for example: --ignore=D101,D2')
        add_check('--convention', metavar='<name>', default=None,
                  help='choose the basic list of checked errors by specifying '
                       'an existing convention. Possible conventions: {}.'
                       .format(', '.join(conventions)))
        add_check('--add-select', metavar='<codes>', default=None,
                  help='add extra error codes to check to the basic list of '
                       'errors previously set by --select, --ignore or '
                       '--convention.')
        add_check('--add-ignore', metavar='<codes>', default=None,
                  help='ignore extra error codes by removing them from the '
                       'basic list previously set by --select, --ignore '
                       'or --convention.')

        parser.add_option_group(check_group)

        # Match clauses
        option('--match', metavar='<pattern>', default=None,
               help=("check only files that exactly match <pattern> regular "
                     "expression; default is --match='{}' which matches "
                     "files that don't start with 'test_' but end with "
                     "'.py'").format(cls.DEFAULT_MATCH_RE))
        option('--match-dir', metavar='<pattern>', default=None,
               help=("search only dirs that exactly match <pattern> regular "
                     "expression; default is --match-dir='{}', which "
                     "matches all dirs that don't start with "
                     "a dot").format(cls.DEFAULT_MATCH_DIR_RE))

        # Decorators
        option('--ignore-decorators', metavar='<decorators>', default=None,
               help=("ignore any functions or methods that are decorated "
                     "by a function with a name fitting the <decorators> "
                     "regular expression; default is --ignore-decorators='{0}'"
                     " which does not ignore any decorated functions."
                     .format(cls.DEFAULT_IGNORE_DECORATORS_RE)))
        return parser


# Check configuration - used by the ConfigurationParser class.
CheckConfiguration = namedtuple('CheckConfiguration',
                                ('checked_codes', 'match', 'match_dir',
                                 'ignore_decorators'))


class IllegalConfiguration(Exception):
    """An exception for illegal configurations."""

    pass


# General configurations for pydocstyle run.
RunConfiguration = namedtuple('RunConfiguration',
                              ('explain', 'source', 'debug',
                               'verbose', 'count', 'config'))
