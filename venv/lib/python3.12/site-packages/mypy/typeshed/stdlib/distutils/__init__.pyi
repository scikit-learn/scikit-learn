# Attempts to improve these stubs are probably not the best use of time:
# - distutils is deleted in Python 3.12 and newer
# - Most users already do not use stdlib distutils, due to setuptools monkeypatching
# - We have very little quality assurance on these stubs, since due to the two above issues
#   we allowlist all distutils errors in stubtest.
