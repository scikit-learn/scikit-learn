import os
import sys
import glob
import importlib

# Utility for testing the presents of pxd-files in the installation
# (and not in ../sklearn-folder)
#
# To add a test for a pxd-file, add an corresponding pyx-file. #
# Attention: See other examples for the naming convertion!
#


# it is ok for build_and_load_pyx to have a deprication warning,
# because pyximport uses imp
def build_and_load_pyx(module_name, include_dirs=[], cpp=False,
                       use_numpy=True, raise_if_pyx_hooks=True):
    """Helper function to build and import a pyx-extension

    Uses pyximport. Deletes PyxImporters from sys.meta_path, if present.

    Parameters
    ----------
    module_name : name of the pyx-module which should be build/imported

    include_dirs : additional include_dirs for C-headers

    cpp : True if should build with c++ rather than c

    use_numpy : True if numpy-headers should be included as well

    raise_if_pyx_hooks : if True raises an error if PyxImporter hooks are
                         already present, rather than deleting them
    """
    import pyximport
    # delete all old hooks or raise if needed:
    old_hooks = [
        hook for hook in sys.meta_path
        if isinstance(hook, pyximport.PyxImporter)
    ]
    if old_hooks and raise_if_pyx_hooks:
        raise RuntimeError('pyx import hooks already installed')
    for pyx_hook in old_hooks:
        pyximport.uninstall(pyx_hook)

    # install new hook:
    script_args = ["--force"]
    if cpp:
        script_args.append("--cython-cplus")
    if use_numpy:
        import numpy as np
        include_dirs.append(np.get_include())
    setup_args = {
        "include_dirs": os.pathsep.join(include_dirs),
        "script_args": script_args,
    }
    py_importer, pyx_importer = pyximport.install(setup_args=setup_args,
                                                  language_level=3)
    try:
        module = importlib.import_module(module_name)
        return module
    finally:
        # clean up the sys.meth_path:
        pyximport.uninstall(py_importer, pyx_importer)


# building and importing all pyx-files for the pattern
def find_and_test(pattern, is_cpp):
    failed = []
    pyx_modules = glob.glob(pattern)
    for file_name in pyx_modules:
        module_name = file_name[0:-4]
        print("\n\ntesting pyx-module", module_name, "...")
        all_ok = False
        try:
            all_ok = build_and_load_pyx(module_name, cpp=is_cpp).all_ok()
        except ImportError:
            pass
        if all_ok:
            print("passed")
        else:
            print("failed")
            failed.append(module_name)
    return failed


# testing c-modules and cpp-modules in this folder
failed = find_and_test('*_c_test.pyx', False)
failed += find_and_test('*_cpp_test.pyx', True)

print("==============================\n\n\n")
if failed:
    print(len(failed), "test(s) failed!")
    print("Failed modules are:", failed)
else:
    print("All ok")
