"""Embed the vcomp dll before generating the scikit-learn Windows wheel.

This script should be run from the root of the scikit-learn source tree,
after running the `python setup.py build` command and before running
the `python setup.py bdist_wheel` command.
"""

import os
import os.path as op
import shutil
from glob import glob
import textwrap

VCOMP140_SRC_PATH = "C:\\Windows\System32\\vcomp140.dll"  # noqa
TARGET_FOLDER_GLOB_PATTERN = "build/lib.*/sklearn"


def make_distributor_init(sklearn_dirname, dll_filename):
    """Create a _distributor_init.py file for the vcomp dll.

    This file is imported first when importing the sklearn package so as
    to pre-load the vendored vcomp dll.
    """
    distributor_init = op.join(sklearn_dirname, '_distributor_init.py')
    with open(distributor_init, 'wt') as f:
        f.write(textwrap.dedent("""
            '''
            Helper to preload the OpenMP dll to prevent "dll not found"
            errors.
            Once a DLL is preloaded, its namespace is made available to any
            subsequent DLL. This file originated in the scikit-learn-wheels
            github repo, and is created as part of the scripts that build the
            wheel.
            '''
            import os
            import os.path as op
            from ctypes import WinDLL


            if os.name == 'nt':
                # Pre-load the DLL stored in sklearn/.libs by convention.
                dll_path = op.join(op.dirname(__file__), '.libs', '{}')
                WinDLL(op.abspath(dll_path))

        """.format(dll_filename)))
    return op.abspath(distributor_init)


def vendor_vcomp140():
    # TODO: use threadpoolctl to dynamically locate the right vcomp dll
    # instead? This would require first in-place building scikit-learn
    # to make it "importable".
    if not op.exists(VCOMP140_SRC_PATH):
        raise ValueError("Could not find %r" % VCOMP140_SRC_PATH)

    if not op.isdir("build"):
        raise RuntimeError("Could not find ./build/ folder. "
                           "Run 'python setup.py build' first")
    target_folders = glob(TARGET_FOLDER_GLOB_PATTERN)
    if len(target_folders) == 0:
        raise RuntimeError("Could not find folder matching '%s'"
                           % TARGET_FOLDER_GLOB_PATTERN)
    if len(target_folders) > 1:
        raise RuntimeError("Found too many target folders: '%s'"
                           % "', '".join(target_folders))
    target_folder = op.abspath(op.join(target_folders[0], ".libs"))

    # create the "sklearn/.libs" subfolder
    if not op.exists(target_folder):
        os.mkdir(target_folder)

    print("Copying '%s' to:\n%s" % (VCOMP140_SRC_PATH, target_folder))
    shutil.copy2(VCOMP140_SRC_PATH, target_folder)

    # Generate the _distributor_init file in the source tree.
    print("Generating the '_distributor_init.py' file in:")
    print(make_distributor_init("sklearn", op.basename(VCOMP140_SRC_PATH)))
