"""Config handling logic for Flake8."""
import collections
import configparser
import logging
import os.path
import sys

from flake8 import utils

LOG = logging.getLogger(__name__)

__all__ = ('ConfigFileFinder', 'MergedConfigParser')


class ConfigFileFinder(object):
    """Encapsulate the logic for finding and reading config files."""

    PROJECT_FILENAMES = ('setup.cfg', 'tox.ini')

    def __init__(self, program_name, args, extra_config_files):
        """Initialize object to find config files.

        :param str program_name:
            Name of the current program (e.g., flake8).
        :param list args:
            The extra arguments passed on the command-line.
        :param list extra_config_files:
            Extra configuration files specified by the user to read.
        """
        # The values of --append-config from the CLI
        extra_config_files = extra_config_files or []
        self.extra_config_files = [
            # Ensure the paths are absolute paths for local_config_files
            os.path.abspath(f) for f in extra_config_files
        ]

        # Platform specific settings
        self.is_windows = sys.platform == 'win32'
        self.xdg_home = os.environ.get('XDG_CONFIG_HOME',
                                       os.path.expanduser('~/.config'))

        # Look for '.<program_name>' files
        self.program_config = '.' + program_name
        self.program_name = program_name

        # List of filenames to find in the local/project directory
        self.project_filenames = ('setup.cfg', 'tox.ini', self.program_config)

        self.local_directory = os.path.abspath(os.curdir)

        if not args:
            args = ['.']
        self.parent = self.tail = os.path.abspath(os.path.commonprefix(args))

        # caches to avoid double-reading config files
        self._local_configs = None
        self._user_config = None
        self._cli_configs = {}

    @staticmethod
    def _read_config(files):
        config = configparser.RawConfigParser()
        if isinstance(files, (str, type(u''))):
            files = [files]

        found_files = []
        for filename in files:
            try:
                found_files.extend(config.read(filename))
            except UnicodeDecodeError:
                LOG.exception("There was an error decoding a config file."
                              "The file with a problem was %s.",
                              filename)
            except configparser.ParsingError:
                LOG.exception("There was an error trying to parse a config "
                              "file. The file with a problem was %s.",
                              filename)
        return (config, found_files)

    def cli_config(self, files):
        """Read and parse the config file specified on the command-line."""
        if files not in self._cli_configs:
            config, found_files = self._read_config(files)
            if found_files:
                LOG.debug('Found cli configuration files: %s', found_files)
            self._cli_configs[files] = config
        return self._cli_configs[files]

    def generate_possible_local_files(self):
        """Find and generate all local config files."""
        tail = self.tail
        parent = self.parent
        found_config_files = False
        while tail and not found_config_files:
            for project_filename in self.project_filenames:
                filename = os.path.abspath(os.path.join(parent,
                                                        project_filename))
                if os.path.exists(filename):
                    yield filename
                    found_config_files = True
                    self.local_directory = parent
            (parent, tail) = os.path.split(parent)

    def local_config_files(self):
        """Find all local config files which actually exist.

        Filter results from
        :meth:`~ConfigFileFinder.generate_possible_local_files` based
        on whether the filename exists or not.

        :returns:
            List of files that exist that are local project config files with
            extra config files appended to that list (which also exist).
        :rtype:
            [str]
        """
        exists = os.path.exists
        return [
            filename
            for filename in self.generate_possible_local_files()
        ] + [f for f in self.extra_config_files if exists(f)]

    def local_configs(self):
        """Parse all local config files into one config object."""
        if self._local_configs is None:
            config, found_files = self._read_config(self.local_config_files())
            if found_files:
                LOG.debug('Found local configuration files: %s', found_files)
            self._local_configs = config
        return self._local_configs

    def user_config_file(self):
        """Find the user-level config file."""
        if self.is_windows:
            return os.path.expanduser('~\\' + self.program_config)
        return os.path.join(self.xdg_home, self.program_name)

    def user_config(self):
        """Parse the user config file into a config object."""
        if self._user_config is None:
            config, found_files = self._read_config(self.user_config_file())
            if found_files:
                LOG.debug('Found user configuration files: %s', found_files)
            self._user_config = config
        return self._user_config


class MergedConfigParser(object):
    """Encapsulate merging different types of configuration files.

    This parses out the options registered that were specified in the
    configuration files, handles extra configuration files, and returns
    dictionaries with the parsed values.
    """

    #: Set of types that should use the
    #: :meth:`~configparser.RawConfigParser.getint` method.
    GETINT_TYPES = {'int', 'count'}
    #: Set of actions that should use the
    #: :meth:`~configparser.RawConfigParser.getbool` method.
    GETBOOL_ACTIONS = {'store_true', 'store_false'}

    def __init__(self, option_manager, config_finder):
        """Initialize the MergedConfigParser instance.

        :param flake8.options.manager.OptionManager option_manager:
            Initialized OptionManager.
        :param flake8.options.config.ConfigFileFinder config_finder:
            Initialized ConfigFileFinder.
        """
        #: Our instance of flake8.options.manager.OptionManager
        self.option_manager = option_manager
        #: The prog value for the cli parser
        self.program_name = option_manager.program_name
        #: Mapping of configuration option names to
        #: :class:`~flake8.options.manager.Option` instances
        self.config_options = option_manager.config_options_dict
        #: Our instance of our :class:`~ConfigFileFinder`
        self.config_finder = config_finder

    def _normalize_value(self, option, value):
        final_value = option.normalize(
            value,
            self.config_finder.local_directory,
        )
        LOG.debug('%r has been normalized to %r for option "%s"',
                  value, final_value, option.config_name)
        return final_value

    def _parse_config(self, config_parser):
        config_dict = {}
        for option_name in config_parser.options(self.program_name):
            if option_name not in self.config_options:
                LOG.debug('Option "%s" is not registered. Ignoring.',
                          option_name)
                continue
            option = self.config_options[option_name]

            # Use the appropriate method to parse the config value
            method = config_parser.get
            if (option.type in self.GETINT_TYPES or
                    option.action in self.GETINT_TYPES):
                method = config_parser.getint
            elif option.action in self.GETBOOL_ACTIONS:
                method = config_parser.getboolean

            value = method(self.program_name, option_name)
            LOG.debug('Option "%s" returned value: %r', option_name, value)

            final_value = self._normalize_value(option, value)
            config_dict[option.config_name] = final_value

        return config_dict

    def is_configured_by(self, config):
        """Check if the specified config parser has an appropriate section."""
        return config.has_section(self.program_name)

    def parse_local_config(self):
        """Parse and return the local configuration files."""
        config = self.config_finder.local_configs()
        if not self.is_configured_by(config):
            LOG.debug('Local configuration files have no %s section',
                      self.program_name)
            return {}

        LOG.debug('Parsing local configuration files.')
        return self._parse_config(config)

    def parse_user_config(self):
        """Parse and return the user configuration files."""
        config = self.config_finder.user_config()
        if not self.is_configured_by(config):
            LOG.debug('User configuration files have no %s section',
                      self.program_name)
            return {}

        LOG.debug('Parsing user configuration files.')
        return self._parse_config(config)

    def parse_cli_config(self, config_path):
        """Parse and return the file specified by --config."""
        config = self.config_finder.cli_config(config_path)
        if not self.is_configured_by(config):
            LOG.debug('CLI configuration files have no %s section',
                      self.program_name)
            return {}

        LOG.debug('Parsing CLI configuration files.')
        return self._parse_config(config)

    def merge_user_and_local_config(self):
        """Merge the parsed user and local configuration files.

        :returns:
            Dictionary of the parsed and merged configuration options.
        :rtype:
            dict
        """
        user_config = self.parse_user_config()
        config = self.parse_local_config()

        for option, value in user_config.items():
            config.setdefault(option, value)

        return config

    def parse(self, cli_config=None, isolated=False):
        """Parse and return the local and user config files.

        First this copies over the parsed local configuration and then
        iterates over the options in the user configuration and sets them if
        they were not set by the local configuration file.

        :param str cli_config:
            Value of --config when specified at the command-line. Overrides
            all other config files.
        :param bool isolated:
            Determines if we should parse configuration files at all or not.
            If running in isolated mode, we ignore all configuration files
        :returns:
            Dictionary of parsed configuration options
        :rtype:
            dict
        """
        if isolated:
            LOG.debug('Refusing to parse configuration files due to user-'
                      'requested isolation')
            return {}

        if cli_config:
            LOG.debug('Ignoring user and locally found configuration files. '
                      'Reading only configuration from "%s" specified via '
                      '--config by the user', cli_config)
            return self.parse_cli_config(cli_config)

        return self.merge_user_and_local_config()


def get_local_plugins(config_finder, cli_config=None, isolated=False):
    """Get local plugins lists from config files.

    :param flake8.options.config.ConfigFileFinder config_finder:
        The config file finder to use.
    :param str cli_config:
        Value of --config when specified at the command-line. Overrides
        all other config files.
    :param bool isolated:
        Determines if we should parse configuration files at all or not.
        If running in isolated mode, we ignore all configuration files
    :returns:
        LocalPlugins namedtuple containing two lists of plugin strings,
        one for extension (checker) plugins and one for report plugins.
    :rtype:
        flake8.options.config.LocalPlugins
    """
    local_plugins = LocalPlugins(extension=[], report=[])
    if isolated:
        LOG.debug('Refusing to look for local plugins in configuration'
                  'files due to user-requested isolation')
        return local_plugins

    if cli_config:
        LOG.debug('Reading local plugins only from "%s" specified via '
                  '--config by the user', cli_config)
        config = config_finder.cli_config(cli_config)
    else:
        config = config_finder.local_configs()

    section = '%s:local-plugins' % config_finder.program_name
    for plugin_type in ['extension', 'report']:
        if config.has_option(section, plugin_type):
            local_plugins_string = config.get(section, plugin_type).strip()
            plugin_type_list = getattr(local_plugins, plugin_type)
            plugin_type_list.extend(utils.parse_comma_separated_list(
                local_plugins_string,
                regexp=utils.LOCAL_PLUGIN_LIST_RE,
            ))
    return local_plugins


LocalPlugins = collections.namedtuple('LocalPlugins', 'extension report')
