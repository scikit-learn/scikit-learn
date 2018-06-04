"""Configuration management setup

Some terminology:
- name
  As written in config files.
- value
  Value associated with a name
- key
  Name combined with it's section (section.name)
- variant
  A single word describing where the configuration key-value pair came from
"""

import locale
import logging
import os

from pip._vendor import six
from pip._vendor.six.moves import configparser

from pip._internal.exceptions import ConfigurationError
from pip._internal.locations import (
    legacy_config_file, new_config_file, running_under_virtualenv,
    site_config_files, venv_config_file,
)
from pip._internal.utils.misc import ensure_dir, enum
from pip._internal.utils.typing import MYPY_CHECK_RUNNING

if MYPY_CHECK_RUNNING:
    from typing import (  # noqa: F401
        Any, Dict, Iterable, List, NewType, Optional, Tuple
    )

    RawConfigParser = configparser.RawConfigParser  # Shorthand
    Kind = NewType("Kind", str)

logger = logging.getLogger(__name__)


# NOTE: Maybe use the optionx attribute to normalize keynames.
def _normalize_name(name):
    # type: (str) -> str
    """Make a name consistent regardless of source (environment or file)
    """
    name = name.lower().replace('_', '-')
    if name.startswith('--'):
        name = name[2:]  # only prefer long opts
    return name


def _disassemble_key(name):
    # type: (str) -> List[str]
    return name.split(".", 1)


# The kinds of configurations there are.
kinds = enum(
    USER="user",        # User Specific
    GLOBAL="global",    # System Wide
    VENV="venv",        # Virtual Environment Specific
    ENV="env",          # from PIP_CONFIG_FILE
    ENV_VAR="env-var",  # from Environment Variables
)


class Configuration(object):
    """Handles management of configuration.

    Provides an interface to accessing and managing configuration files.

    This class converts provides an API that takes "section.key-name" style
    keys and stores the value associated with it as "key-name" under the
    section "section".

    This allows for a clean interface wherein the both the section and the
    key-name are preserved in an easy to manage form in the configuration files
    and the data stored is also nice.
    """

    def __init__(self, isolated, load_only=None):
        # type: (bool, Kind) -> None
        super(Configuration, self).__init__()

        _valid_load_only = [kinds.USER, kinds.GLOBAL, kinds.VENV, None]
        if load_only not in _valid_load_only:
            raise ConfigurationError(
                "Got invalid value for load_only - should be one of {}".format(
                    ", ".join(map(repr, _valid_load_only[:-1]))
                )
            )
        self.isolated = isolated  # type: bool
        self.load_only = load_only  # type: Optional[Kind]

        # The order here determines the override order.
        self._override_order = [
            kinds.GLOBAL, kinds.USER, kinds.VENV, kinds.ENV, kinds.ENV_VAR
        ]

        self._ignore_env_names = ["version", "help"]

        # Because we keep track of where we got the data from
        self._parsers = {
            variant: [] for variant in self._override_order
        }  # type: Dict[Kind, List[Tuple[str, RawConfigParser]]]
        self._config = {
            variant: {} for variant in self._override_order
        }  # type: Dict[Kind, Dict[str, Any]]
        self._modified_parsers = []  # type: List[Tuple[str, RawConfigParser]]

    def load(self):
        # type: () -> None
        """Loads configuration from configuration files and environment
        """
        self._load_config_files()
        if not self.isolated:
            self._load_environment_vars()

    def get_file_to_edit(self):
        # type: () -> Optional[str]
        """Returns the file with highest priority in configuration
        """
        assert self.load_only is not None, \
            "Need to be specified a file to be editing"

        try:
            return self._get_parser_to_modify()[0]
        except IndexError:
            return None

    def items(self):
        # type: () -> Iterable[Tuple[str, Any]]
        """Returns key-value pairs like dict.items() representing the loaded
        configuration
        """
        return self._dictionary.items()

    def get_value(self, key):
        # type: (str) -> Any
        """Get a value from the configuration.
        """
        try:
            return self._dictionary[key]
        except KeyError:
            raise ConfigurationError("No such key - {}".format(key))

    def set_value(self, key, value):
        # type: (str, Any) -> None
        """Modify a value in the configuration.
        """
        self._ensure_have_load_only()

        fname, parser = self._get_parser_to_modify()

        if parser is not None:
            section, name = _disassemble_key(key)

            # Modify the parser and the configuration
            if not parser.has_section(section):
                parser.add_section(section)
            parser.set(section, name, value)

        self._config[self.load_only][key] = value
        self._mark_as_modified(fname, parser)

    def unset_value(self, key):
        # type: (str) -> None
        """Unset a value in the configuration.
        """
        self._ensure_have_load_only()

        if key not in self._config[self.load_only]:
            raise ConfigurationError("No such key - {}".format(key))

        fname, parser = self._get_parser_to_modify()

        if parser is not None:
            section, name = _disassemble_key(key)

            # Remove the key in the parser
            modified_something = False
            if parser.has_section(section):
                # Returns whether the option was removed or not
                modified_something = parser.remove_option(section, name)

            if modified_something:
                # name removed from parser, section may now be empty
                section_iter = iter(parser.items(section))
                try:
                    val = six.next(section_iter)
                except StopIteration:
                    val = None

                if val is None:
                    parser.remove_section(section)

                self._mark_as_modified(fname, parser)
            else:
                raise ConfigurationError(
                    "Fatal Internal error [id=1]. Please report as a bug."
                )

        del self._config[self.load_only][key]

    def save(self):
        # type: () -> None
        """Save the currentin-memory state.
        """
        self._ensure_have_load_only()

        for fname, parser in self._modified_parsers:
            logger.info("Writing to %s", fname)

            # Ensure directory exists.
            ensure_dir(os.path.dirname(fname))

            with open(fname, "w") as f:
                parser.write(f)  # type: ignore

    #
    # Private routines
    #

    def _ensure_have_load_only(self):
        # type: () -> None
        if self.load_only is None:
            raise ConfigurationError("Needed a specific file to be modifying.")
        logger.debug("Will be working with %s variant only", self.load_only)

    @property
    def _dictionary(self):
        # type: () -> Dict[str, Any]
        """A dictionary representing the loaded configuration.
        """
        # NOTE: Dictionaries are not populated if not loaded. So, conditionals
        #       are not needed here.
        retval = {}

        for variant in self._override_order:
            retval.update(self._config[variant])

        return retval

    def _load_config_files(self):
        # type: () -> None
        """Loads configuration from configuration files
        """
        config_files = dict(self._iter_config_files())
        if config_files[kinds.ENV][0:1] == [os.devnull]:
            logger.debug(
                "Skipping loading configuration files due to "
                "environment's PIP_CONFIG_FILE being os.devnull"
            )
            return

        for variant, files in config_files.items():
            for fname in files:
                # If there's specific variant set in `load_only`, load only
                # that variant, not the others.
                if self.load_only is not None and variant != self.load_only:
                    logger.debug(
                        "Skipping file '%s' (variant: %s)", fname, variant
                    )
                    continue

                parser = self._load_file(variant, fname)

                # Keeping track of the parsers used
                self._parsers[variant].append((fname, parser))

    def _load_file(self, variant, fname):
        # type: (Kind, str) -> RawConfigParser
        logger.debug("For variant '%s', will try loading '%s'", variant, fname)
        parser = self._construct_parser(fname)

        for section in parser.sections():
            items = parser.items(section)
            self._config[variant].update(self._normalized_keys(section, items))

        return parser

    def _construct_parser(self, fname):
        # type: (str) -> RawConfigParser
        parser = configparser.RawConfigParser()
        # If there is no such file, don't bother reading it but create the
        # parser anyway, to hold the data.
        # Doing this is useful when modifying and saving files, where we don't
        # need to construct a parser.
        if os.path.exists(fname):
            try:
                parser.read(fname)
            except UnicodeDecodeError:
                raise ConfigurationError((
                    "ERROR: "
                    "Configuration file contains invalid %s characters.\n"
                    "Please fix your configuration, located at %s\n"
                ) % (locale.getpreferredencoding(False), fname))
        return parser

    def _load_environment_vars(self):
        # type: () -> None
        """Loads configuration from environment variables
        """
        self._config[kinds.ENV_VAR].update(
            self._normalized_keys(":env:", self._get_environ_vars())
        )

    def _normalized_keys(self, section, items):
        # type: (str, Iterable[Tuple[str, Any]]) -> Dict[str, Any]
        """Normalizes items to construct a dictionary with normalized keys.

        This routine is where the names become keys and are made the same
        regardless of source - configuration files or environment.
        """
        normalized = {}
        for name, val in items:
            key = section + "." + _normalize_name(name)
            normalized[key] = val
        return normalized

    def _get_environ_vars(self):
        # type: () -> Iterable[Tuple[str, str]]
        """Returns a generator with all environmental vars with prefix PIP_"""
        for key, val in os.environ.items():
            should_be_yielded = (
                key.startswith("PIP_") and
                key[4:].lower() not in self._ignore_env_names
            )
            if should_be_yielded:
                yield key[4:].lower(), val

    # XXX: This is patched in the tests.
    def _iter_config_files(self):
        # type: () -> Iterable[Tuple[Kind, List[str]]]
        """Yields variant and configuration files associated with it.

        This should be treated like items of a dictionary.
        """
        # SMELL: Move the conditions out of this function

        # environment variables have the lowest priority
        config_file = os.environ.get('PIP_CONFIG_FILE', None)
        if config_file is not None:
            yield kinds.ENV, [config_file]
        else:
            yield kinds.ENV, []

        # at the base we have any global configuration
        yield kinds.GLOBAL, list(site_config_files)

        # per-user configuration next
        should_load_user_config = not self.isolated and not (
            config_file and os.path.exists(config_file)
        )
        if should_load_user_config:
            # The legacy config file is overridden by the new config file
            yield kinds.USER, [legacy_config_file, new_config_file]

        # finally virtualenv configuration first trumping others
        if running_under_virtualenv():
            yield kinds.VENV, [venv_config_file]

    def _get_parser_to_modify(self):
        # type: () -> Tuple[str, RawConfigParser]
        # Determine which parser to modify
        parsers = self._parsers[self.load_only]
        if not parsers:
            # This should not happen if everything works correctly.
            raise ConfigurationError(
                "Fatal Internal error [id=2]. Please report as a bug."
            )

        # Use the highest priority parser.
        return parsers[-1]

    # XXX: This is patched in the tests.
    def _mark_as_modified(self, fname, parser):
        # type: (str, RawConfigParser) -> None
        file_parser_tuple = (fname, parser)
        if file_parser_tuple not in self._modified_parsers:
            self._modified_parsers.append(file_parser_tuple)
