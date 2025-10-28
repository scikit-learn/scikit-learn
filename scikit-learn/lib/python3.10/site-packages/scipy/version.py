
"""
Module to expose more detailed version info for the installed `scipy`
"""
version = "1.15.3"
full_version = version
short_version = version.split('.dev')[0]
git_revision = "e29dcb65a2040f04819b426a04b60d44a8f69c04"
release = 'dev' not in version and '+' not in version

if not release:
    version = full_version
