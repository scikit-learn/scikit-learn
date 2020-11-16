"""Embed VCRUNTIME140.dll after generating the Windows wheels."""


import os
import os.path as op
import shutil
import sys
import textwrap


TARGET_FOLDER = op.join("sklearn", ".libs")
DISTRIBUTOR_INIT = op.join("sklearn", "_distributor_init.py")
VCRUNTIME140_SRC_PATH = "C:\\Windows\\System32\\VCRUNTIME140.dll"


def make_distributor_init(distributor_init, vcruntime140):
    """Create a _distributor_init.py file for VCRUNTIME140.dll.

    This file is imported first when importing the sklearn
    package so as to pre-load the vendored VCRUNTIME140.dll.
    """
    with open(distributor_init, "wt") as f:
        f.write(textwrap.dedent("""
            '''Helper to preload VCRUNTIME140.dll to prevent "not found"
            errors.

            Once VCRUNTIME140.dll is preloaded, the namespace is made
            available to any subsequent VCRUNTIME140.dll. This is
            created as part of the scripts that build the wheel.
            '''


            import os
            import os.path as op
            from ctypes import WinDLL


            if os.name == "nt":
                # Load VCRUNTIME140.dll in sklearn/.libs
                vcruntime140 = op.join(op.dirname(__file__), ".libs", "{0}")
                WinDLL(op.abspath(vcruntime140))
            """.format(vcruntime140)))


def main(wheel_dirname):
    """Embed VCRUNTIME140.dll in the wheel."""
    if not op.exists(VCRUNTIME140_SRC_PATH):
        raise ValueError(f"Could not find {VCRUNTIME140_SRC_PATH}.")

    if not op.isdir(wheel_dirname):
        raise RuntimeError(f"Could not find {wheel_dirname} file.")

    vcruntime140_dll = op.basename(VCRUNTIME140_SRC_PATH)
    target_folder = op.join(wheel_dirname, TARGET_FOLDER)
    distributor_init = op.join(wheel_dirname, DISTRIBUTOR_INIT)

    # Create the "sklearn/.libs" subfolder
    if not op.exists(target_folder):
        os.mkdir(target_folder)

    print(f"Copying {VCRUNTIME140_SRC_PATH} to {target_folder}.")
    shutil.copy2(VCRUNTIME140_SRC_PATH, target_folder)

    # Generate the _distributor_init file in the source tree
    print("Generating the '_distributor_init.py' file.")
    make_distributor_init(distributor_init, vcruntime140_dll)


if __name__ == "__main__":
    _, wheel_file = sys.argv
    main(wheel_file)
