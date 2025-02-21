from __future__ import absolute_import
import sys

if sys.version_info < (3, 5):
    # _pyximport3 module requires at least Python 3.5
    from pyximport._pyximport2 import install, uninstall, show_docs
else:
    from pyximport._pyximport3 import install, uninstall, show_docs

if __name__ == '__main__':
    show_docs()
