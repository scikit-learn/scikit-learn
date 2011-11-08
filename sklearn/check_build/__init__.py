import os

try:
    from ._check_build import check_build
except ImportError, e:
    # Raise a comprehensible error and list the contents of the
    # directory to help debugging on the mailing list.
    local_dir = os.path.split(__file__)[0]
    dir_content = list()
    for i, filename in enumerate(os.listdir(local_dir)):
        if ((i + 1) % 3):
            dir_content.append(filename.ljust(26))
        else:
            dir_content.append(filename + '\n')
    raise ImportError(
"""%s
___________________________________________________________________________
Contents of %s:
%s
___________________________________________________________________________
It seems that the scikit-learn has not been built correctly.

If you have installed the scikit-learn from source, please do not forget
to build the package before using it: run `python setup.py install` or
`make` in the source directory.

If you have used an installer, please check that it is suited for your
Python version, your operating system and your platform.
""" % (e, local_dir, ''.join(dir_content).strip()))


