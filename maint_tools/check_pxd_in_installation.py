"""Utility for testing presence and usability of .pxd files in the installation

Usage:
------
python check_pxd_in_installation.py path/to/install_dir/of/scikit-learn
"""

import os
import sys
import glob
import tempfile
import textwrap
import subprocess


skl_dir = sys.argv[1]
skl_dir = os.path.abspath(skl_dir)

# find all .pxd files
pxd_files = glob.glob(os.path.join(skl_dir, '**', '*.pxd'), recursive=True)
pxd_files = [f.replace(skl_dir, 'sklearn') for f in pxd_files]

print("> Found pxd files:")
for file in pxd_files:
    print(' -', file)

print("\n> Trying to compile a cython extension cimporting all corresponding "
      "modules\n")
start_dir = os.path.abspath('.')

with tempfile.TemporaryDirectory() as tmpdir:
    try:
        os.chdir(tmpdir)

        # A cython test file which cimports all modules corresponding to found
        # pxd files.
        # e.g. sklearn/tree/_utils.pxd becomes `cimport sklearn.tree._utils`
        with open('tst.pyx', 'w') as f:
            for file in pxd_files:
                to_import = file.replace(os.path.sep, '.')
                to_import = to_import.replace('.pxd', '')
                f.write('cimport ' + to_import + '\n')

        # A basic setup file to build the test file.
        # We set the language to c++ and we use numpy.get_include() because
        # some modules require it.
        with open('setup_tst.py', 'w') as f:
            f.write(textwrap.dedent(
                """
                from distutils.core import setup
                from distutils.extension import Extension
                from Cython.Build import cythonize
                import numpy

                extensions = [Extension("tst",
                                        sources=["tst.pyx"],
                                        language="c++",
                                        include_dirs=[numpy.get_include()])]

                setup(ext_modules=cythonize(extensions))
                """))

        subprocess.run(["python", "setup_tst.py", "build_ext", "-i"],
                       check=True)

    except subprocess.CalledProcessError:
        raise

    finally:
        os.chdir(start_dir)

    print("\n> Compilation succeeded !")
