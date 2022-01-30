import sys

if sys.platform == "win32":
    from ntpath import *
else:
    from posixpath import *
