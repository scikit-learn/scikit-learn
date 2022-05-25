"""Embed vcomp140.dll, vcruntime140.dll and vcruntime140_1.dll.

Note that vcruntime140_1.dll is only required (and available)
for 64-bit architectures.
"""


import os
import os.path as op
import shutil
import sys
import textwrap


TARGET_FOLDER = op.join("sklearn", ".libs")
DISTRIBUTOR_INIT = op.join("sklearn", "_distributor_init.py")
VCOMP140_SRC_PATH = "C:\\Windows\\System32\\vcomp140.dll"
VCRUNTIME140_SRC_PATH = "C:\\Windows\\System32\\vcruntime140.dll"
VCRUNTIME140_1_SRC_PATH = "C:\\Windows\\System32\\vcruntime140_1.dll"


def make_distributor_init_32_bits(
    distributor_init, vcomp140_dll_filename, vcruntime140_dll_filename
):
    """Create a _distributor_init.py file for 32-bit architectures.

    This file is imported first when importing the sklearn package
    so as to pre-load the vendored vcomp140.dll and vcruntime140.dll.
    """
    with open(distributor_init, "wt") as f:
        f.write(
            textwrap.dedent(
                """
            '''Helper to preload vcomp140.dll and vcruntime140.dll to
            prevent "not found" errors.

            Once vcomp140.dll and vcruntime140.dll are preloaded, the
            namespace is made available to any subsequent vcomp140.dll
            and vcruntime140.dll. This is created as part of the scripts
            that build the wheel.
            '''


            import os
            import os.path as op
            from ctypes import WinDLL


            if os.name == "nt":
                # Load vcomp140.dll and vcruntime140.dll
                libs_path = op.join(op.dirname(__file__), ".libs")
                vcomp140_dll_filename = op.join(libs_path, "{0}")
                vcruntime140_dll_filename = op.join(libs_path, "{1}")
                WinDLL(op.abspath(vcomp140_dll_filename))
                WinDLL(op.abspath(vcruntime140_dll_filename))
            """.format(
                    vcomp140_dll_filename, vcruntime140_dll_filename
                )
            )
        )


def make_distributor_init_64_bits(
    distributor_init,
    vcomp140_dll_filename,
    vcruntime140_dll_filename,
    vcruntime140_1_dll_filename,
):
    """Create a _distributor_init.py file for 64-bit architectures.

    This file is imported first when importing the sklearn package
    so as to pre-load the vendored vcomp140.dll, vcruntime140.dll
    and vcruntime140_1.dll.
    """
    with open(distributor_init, "wt") as f:
        f.write(
            textwrap.dedent(
                """
            '''Helper to preload vcomp140.dll, vcruntime140.dll and
            vcruntime140_1.dll to prevent "not found" errors.

            Once vcomp140.dll, vcruntime140.dll and vcruntime140_1.dll are
            preloaded, the namespace is made available to any subsequent
            vcomp140.dll, vcruntime140.dll and vcruntime140_1.dll. This is
            created as part of the scripts that build the wheel.
            '''


            import os
            import os.path as op
            from ctypes import WinDLL


            if os.name == "nt":
                # Load vcomp140.dll, vcruntime140.dll and vcruntime140_1.dll
                libs_path = op.join(op.dirname(__file__), ".libs")
                vcomp140_dll_filename = op.join(libs_path, "{0}")
                vcruntime140_dll_filename = op.join(libs_path, "{1}")
                vcruntime140_1_dll_filename = op.join(libs_path, "{2}")
                WinDLL(op.abspath(vcomp140_dll_filename))
                WinDLL(op.abspath(vcruntime140_dll_filename))
                WinDLL(op.abspath(vcruntime140_1_dll_filename))
            """.format(
                    vcomp140_dll_filename,
                    vcruntime140_dll_filename,
                    vcruntime140_1_dll_filename,
                )
            )
        )


def main(wheel_dirname, bitness):
    """Embed vcomp140.dll, vcruntime140.dll and vcruntime140_1.dll."""
    if not op.exists(VCOMP140_SRC_PATH):
        raise ValueError(f"Could not find {VCOMP140_SRC_PATH}.")

    if not op.exists(VCRUNTIME140_SRC_PATH):
        raise ValueError(f"Could not find {VCRUNTIME140_SRC_PATH}.")

    if not op.exists(VCRUNTIME140_1_SRC_PATH) and bitness == "64":
        raise ValueError(f"Could not find {VCRUNTIME140_1_SRC_PATH}.")

    if not op.isdir(wheel_dirname):
        raise RuntimeError(f"Could not find {wheel_dirname} file.")

    vcomp140_dll_filename = op.basename(VCOMP140_SRC_PATH)
    vcruntime140_dll_filename = op.basename(VCRUNTIME140_SRC_PATH)
    vcruntime140_1_dll_filename = op.basename(VCRUNTIME140_1_SRC_PATH)

    target_folder = op.join(wheel_dirname, TARGET_FOLDER)
    distributor_init = op.join(wheel_dirname, DISTRIBUTOR_INIT)

    # Create the "sklearn/.libs" subfolder
    if not op.exists(target_folder):
        os.mkdir(target_folder)

    print(f"Copying {VCOMP140_SRC_PATH} to {target_folder}.")
    shutil.copy2(VCOMP140_SRC_PATH, target_folder)

    print(f"Copying {VCRUNTIME140_SRC_PATH} to {target_folder}.")
    shutil.copy2(VCRUNTIME140_SRC_PATH, target_folder)

    if bitness == "64":
        print(f"Copying {VCRUNTIME140_1_SRC_PATH} to {target_folder}.")
        shutil.copy2(VCRUNTIME140_1_SRC_PATH, target_folder)

    # Generate the _distributor_init file in the source tree
    print("Generating the '_distributor_init.py' file.")
    if bitness == "32":
        make_distributor_init_32_bits(
            distributor_init, vcomp140_dll_filename, vcruntime140_dll_filename
        )
    else:
        make_distributor_init_64_bits(
            distributor_init,
            vcomp140_dll_filename,
            vcruntime140_dll_filename,
            vcruntime140_1_dll_filename,
        )


if __name__ == "__main__":
    _, wheel_file, bitness = sys.argv
    main(wheel_file, bitness)
