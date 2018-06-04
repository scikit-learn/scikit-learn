# Copyright 2017 Palantir Technologies, Inc.
import configparser
import logging
import os
import sys

log = logging.getLogger(__name__)


class ConfigSource(object):
    """Base class for implementing a config source."""

    def __init__(self, root_path):
        self.root_path = root_path
        self.is_windows = sys.platform == 'win32'
        self.xdg_home = os.environ.get(
            'XDG_CONFIG_HOME', os.path.expanduser('~/.config')
        )

        self._modified_times = {}
        self._configs_cache = {}

    def user_config(self):
        """Return user-level (i.e. home directory) configuration."""
        raise NotImplementedError()

    def project_config(self, document_path):
        """Return project-level (i.e. workspace directory) configuration."""
        raise NotImplementedError()

    def read_config_from_files(self, files):
        files = tuple([f for f in files if os.path.exists(f) and not os.path.isdir(f)])
        modified = tuple([os.path.getmtime(f) for f in files])

        if files in self._modified_times and modified == self._modified_times[files]:
            log.debug("Using cached configuration for %s", files)
            return self._configs_cache[files]

        config = configparser.RawConfigParser()
        found_files = []
        for filename in files:
            found_files.extend(config.read(filename))

        self._configs_cache[files] = config
        self._modified_times[files] = modified

        return config

    @staticmethod
    def parse_config(config, key, options):
        """Parse the config with the given options."""
        conf = {}
        for source, destination, opt_type in options:
            opt_value = _get_opt(config, key, source, opt_type)
            if opt_value is not None:
                _set_opt(conf, destination, opt_value)
        return conf


def _get_opt(config, key, option, opt_type):
    """Get an option from a configparser with the given type."""
    for opt_key in [option, option.replace('-', '_')]:
        if not config.has_option(key, opt_key):
            continue

        if opt_type == bool:
            return config.getbool(key, opt_key)

        if opt_type == int:
            return config.getint(key, opt_key)

        if opt_type == str:
            return config.get(key, opt_key)

        if opt_type == list:
            return _parse_list_opt(config.get(key, opt_key))

        raise ValueError("Unknown option type: %s" % opt_type)


def _parse_list_opt(string):
    return [s.strip() for s in string.split(",") if s.strip()]


def _set_opt(config_dict, path, value):
    """Set the value in the dictionary at the given path if the value is not None."""
    if value is None:
        return

    if '.' not in path:
        config_dict[path] = value
        return

    key, rest = path.split(".", 1)
    if key not in config_dict:
        config_dict[key] = {}

    _set_opt(config_dict[key], rest, value)
