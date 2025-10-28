"""Subpackage to handle settings parsing."""

from __future__ import annotations

from towncrier._settings import load


load_config = load.load_config
ConfigError = load.ConfigError
load_config_from_options = load.load_config_from_options

# Help message for --config CLI option, shared by all sub-commands.
config_option_help = (
    "Pass a custom config file at FILE_PATH. "
    "Default: towncrier.toml or pyproject.toml file, "
    "if both files exist, the first will take precedence."
)

__all__ = [
    "config_option_help",
    "load_config",
    "ConfigError",
    "load_config_from_options",
]
