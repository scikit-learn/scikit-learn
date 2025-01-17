
"""
Module to expose more detailed version info for the installed `numpy`
"""
version = "2.2.1"
__version__ = version
full_version = version

git_revision = "7469245786b57405b7ae264b386400061b3b25d3"
release = 'dev' not in version and '+' not in version
short_version = version.split("+")[0]
