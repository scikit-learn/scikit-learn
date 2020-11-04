"""Embed vcomp140.dll after generating the scikit-learn Windows wheel."""


import os
import os.path as op
import zipfile
import textwrap
from glob import glob


VCOMP140_SRC_PATH = "C:\\Windows\System32\\vcomp140.dll"  # noqa
GLOB_PATTERN = "*.whl"


def make_distributor_init(zipf, sklearn_dirname, dll_filename):
    """Create a _distributor_init.py file for the vcomp140.dll.

    This file is imported first when importing the
    sklearn package so as to pre-load the vendored
    vcomp140.dll.
    """
    distributor_init = op.join(sklearn_dirname, "_distributor_init.py")

    # Remove the _distributor_init file to avoid duplicated files
    os.system(f"zip -d {zipf.filename} \"{distributor_init}\"")

    zipf.writestr(distributor_init, textwrap.dedent("""
        '''Helper to preload vcomp140.dll to prevent "not found" errors.

        Once the vcomp140.dll is preloaded, the namespace is made
        available to any subsequent vcomp140.dll. This is created
        as part of the scripts that build the wheel.
        '''

        import os
        import os.path as op
        from ctypes import WinDLL


        if os.name == "nt":
            # Load the vcomp140.dll in sklearn/.libs by convention
            dll_path = op.join(op.dirname(__file__), ".libs", "{0}")
            WinDLL(op.abspath(dll_path))
        """.format(dll_filename)))


def embed_vcomp140(wheel_dirname):
    """Embed the vcomp140.dll in the wheel."""
    if not op.exists(VCOMP140_SRC_PATH):
        raise ValueError(f"Could not find {VCOMP140_SRC_PATH}.")

    if not op.isdir(wheel_dirname):
        raise RuntimeError(f"Could not find {wheel_dirname} folder.")

    wheel_file_glob_pattern = op.join(wheel_dirname, GLOB_PATTERN)
    wheel_files = glob(wheel_file_glob_pattern)

    if len(wheel_file_glob_pattern) == 0:
        raise RuntimeError(f"Could not find wheels matching "
                           f"{wheel_file_glob_pattern}.")

    if len(wheel_files) > 1:
        raise RuntimeError(f"Found too many wheels: "
                           f"{', '.join(wheel_files)}.")

    wheel_file = wheel_files[0]
    dll_filename = op.basename(VCOMP140_SRC_PATH)
    target_folder = op.join("sklearn", ".libs", "vcomp140.dll")

    with zipfile.ZipFile(wheel_file, "a") as zipf:
        print(f"Copying {VCOMP140_SRC_PATH} to {target_folder}.")
        zipf.write(VCOMP140_SRC_PATH, target_folder)

        # Generate the _distributor_init file in the source tree
        print("Generating the '_distributor_init.py' file.")
        make_distributor_init(zipf, "sklearn", dll_filename)
