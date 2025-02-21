
"""
Module to expose more detailed version info for the installed `scipy`
"""
version = "1.15.1"
full_version = version
short_version = version.split('.dev')[0]
git_revision = "df134eab5a500c2146ed4552c8674a78d8154ee9"
release = 'dev' not in version and '+' not in version

if not release:
    version = full_version
