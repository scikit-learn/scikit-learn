"""A lil' TOML parser."""

__all__ = ("loads", "load", "TOMLDecodeError")
__version__ = "2.0.0"  # DO NOT EDIT THIS LINE MANUALLY. LET bump2version UTILITY DO IT

from tomli._parser import TOMLDecodeError, load, loads

# Pretend this exception was created here.
TOMLDecodeError.__module__ = __name__
