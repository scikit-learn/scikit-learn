""" Distributor init file

Distributors: you can add custom code here to support particular distributions
of scipy.

For example, this is a good place to put any checks for hardware requirements.

The scipy standard source distribution will not put code in this file, so you
can safely replace this file with your own version.
"""

import os

# on Windows SciPy loads important DLLs
# and the code below aims to alleviate issues with DLL
# path resolution portability with an absolute path DLL load
if os.name == 'nt':
    from ctypes import WinDLL
    import glob
    # convention for storing / loading the DLL from
    # scipy/.libs/, if present
    libs_path = os.path.abspath(os.path.join(os.path.dirname(__file__),
                                             '.libs'))
    if os.path.isdir(libs_path):
        # for Python >= 3.8, DLL resolution ignores the PATH variable
        # and the current working directory; see release notes:
        # https://docs.python.org/3/whatsnew/3.8.html#bpo-36085-whatsnew
        # Only the system paths, the directory containing the DLL, and 
        # directories added with add_dll_directory() are searched for 
        # load-time dependencies with Python >= 3.8

        # this module was originally added to support DLL resolution in
        # Python 3.8 because of the changes described above--providing the
        # absolute paths to the DLLs allowed for proper resolution/loading

        # however, we also started to receive reports of problems with DLL
        # resolution with Python 3.7 that were sometimes alleviated with
        # inclusion of the _distributor_init.py module; see SciPy main 
        # repo gh-11826

        # we noticed in scipy-wheels repo gh-70 that inclusion of
        # _distributor_init.py in 32-bit wheels for Python 3.7 resulted
        # in failures in DLL resolution (64-bit 3.7 did not)
        # as a result, we decided to combine both the old (working directory)
        # and new (absolute path to DLL location) DLL resolution mechanisms
        # to improve the chances of resolving DLLs across a wider range of
        # Python versions

        # we did not experiment with manipulating the PATH environment variable
        # to include libs_path; it is not immediately clear if this would have
        # robustness or security advantages over changing working directories
        # as done below

        # we should remove the working directory shims when our minimum supported
        # Python version is 3.8 and trust the improvements to secure DLL loading
        # in the standard lib for Python >= 3.8
        try:
            owd = os.getcwd()
            os.chdir(libs_path)
            for filename in glob.glob(os.path.join(libs_path, '*dll')):
                WinDLL(os.path.abspath(filename))
        finally:
            os.chdir(owd)
