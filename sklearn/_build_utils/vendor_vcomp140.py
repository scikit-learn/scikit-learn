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
GLOB_PATTERN = "lib.*/sklearn"


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


def embed_vcomp140(build_dirname):
    # TODO: use threadpoolctl to dynamically locate the right vcomp dll
    # instead? This would require first in-place building scikit-learn
    # to make it "importable".
    if not op.exists(VCOMP140_SRC_PATH):
        raise ValueError(f"Could not find {VCOMP140_SRC_PATH}.")

    if not op.isdir(build_dirname):
        raise RuntimeError(f"Could not find {build_dirname} folder. "
                           "Run 'python setup.py build' first.")

    target_folder_glob_pattern = op.join(build_dirname, GLOB_PATTERN)
    target_folders = glob(target_folder_glob_pattern)

    if len(target_folders) == 0:
        raise RuntimeError(f"Could not find folder matching "
                           f"{target_folder_glob_pattern}.")

    if len(target_folders) > 1:
        raise RuntimeError(f"Found too many target folders: "
                           f"{', '.join(target_folders)}.")

    target_folder = op.abspath(op.join(target_folders[0], ".libs"))

    # create the "sklearn/.libs" subfolder
    if not op.exists(target_folder):
        os.mkdir(target_folder)

    print("Copying '%s' to:\n%s" % (VCOMP140_SRC_PATH, target_folder))
    shutil.copy2(VCOMP140_SRC_PATH, target_folder)

    # Generate the _distributor_init file in the source tree.
    print("Generating the '_distributor_init.py' file in:")
    print(make_distributor_init("sklearn", op.basename(VCOMP140_SRC_PATH)))
