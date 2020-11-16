"""Embed vcomp140.dll after generating the scikit-learn Windows wheel."""


import os
import os.path as op
import shutil
import sys
import textwrap


TARGET_FOLDER = op.join("sklearn", ".libs")
DISTRIBUTOR_INIT = op.join("sklearn", "_distributor_init.py")
VCOMP140_SRC_PATH = "C:\\Windows\System32\\vcomp140.dll"  # noqa


def make_distributor_init(distributor_init, dll_filename):
    """Create a _distributor_init.py file for the vcomp140.dll.

    This file is imported first when importing the
    sklearn package so as to pre-load the vendored
    vcomp140.dll.
    """
    with open(distributor_init, "wt") as f:
        f.write(textwrap.dedent("""
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


def main(wheel_dirname):
    """Embed the vcomp140.dll in the wheel."""
    if not op.exists(VCOMP140_SRC_PATH):
        raise ValueError(f"Could not find {VCOMP140_SRC_PATH}.")

    if not op.isdir(wheel_dirname):
        raise RuntimeError(f"Could not find {wheel_dirname} file.")

    dll_filename = op.basename(VCOMP140_SRC_PATH)
    target_folder = op.join(wheel_dirname, TARGET_FOLDER)
    distributor_init = op.join(wheel_dirname, DISTRIBUTOR_INIT)

    # Create the "sklearn/.libs" subfolder
    if not op.exists(target_folder):
        os.mkdir(target_folder)

    print(f"Copying {VCOMP140_SRC_PATH} to {target_folder}.")
    shutil.copy2(VCOMP140_SRC_PATH, target_folder)

    # Generate the _distributor_init file in the source tree
    print("Generating the '_distributor_init.py' file.")
    make_distributor_init(distributor_init, dll_filename)


if __name__ == "__main__":
    _, wheel_file = sys.argv
    main(wheel_file)
