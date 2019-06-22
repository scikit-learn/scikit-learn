import re
import os
import glob
from distutils.dep_util import newer


__all__ = ['needs_g77_abi_wrapper', 'split_fortran_files',
           'get_g77_abi_wrappers']


def uses_mkl(info):
    r_mkl = re.compile("mkl")
    libraries = info.get('libraries', '')
    for library in libraries:
        if r_mkl.search(library):
            return True

    return False


def needs_g77_abi_wrapper(info):
    """Returns True if g77 ABI wrapper must be used."""
    return uses_mkl(info)


def get_g77_abi_wrappers(info):
    """
    Returns file names of source files containing Fortran ABI wrapper
    routines.
    """
    wrapper_sources = []

    path = os.path.abspath(os.path.dirname(__file__))
    if needs_g77_abi_wrapper(info):
        wrapper_sources += [
            os.path.join(path, 'src', 'wrap_g77_abi_f.f'),
            os.path.join(path, 'src', 'wrap_g77_abi_c.c'),
        ]
    else:
        wrapper_sources += [
            os.path.join(path, 'src', 'wrap_dummy_g77_abi.f'),
        ]
    return wrapper_sources


def split_fortran_files(source_dir, subroutines=None):
    """Split each file in `source_dir` into separate files per subroutine.

    Parameters
    ----------
    source_dir : str
        Full path to directory in which sources to be split are located.
    subroutines : list of str, optional
        Subroutines to split. (Default: all)

    Returns
    -------
    fnames : list of str
        List of file names (not including any path) that were created
        in `source_dir`.

    Notes
    -----
    This function is useful for code that can't be compiled with g77 because of
    type casting errors which do work with gfortran.

    Created files are named: ``original_name + '_subr_i' + '.f'``, with ``i``
    starting at zero and ending at ``num_subroutines_in_file - 1``.

    """

    if subroutines is not None:
        subroutines = [x.lower() for x in subroutines]

    def split_file(fname):
        with open(fname, 'rb') as f:
            lines = f.readlines()
            subs = []
            need_split_next = True

            # find lines with SUBROUTINE statements
            for ix, line in enumerate(lines):
                m = re.match(b'^\\s+subroutine\\s+([a-z0-9_]+)\\s*\\(', line, re.I)
                if m and line[0] not in b'Cc!*':
                    if subroutines is not None:
                        subr_name = m.group(1).decode('ascii').lower()
                        subr_wanted = (subr_name in subroutines)
                    else:
                        subr_wanted = True
                    if subr_wanted or need_split_next:
                        need_split_next = subr_wanted
                        subs.append(ix)

            # check if no split needed
            if len(subs) <= 1:
                return [fname]

            # write out one file per subroutine
            new_fnames = []
            num_files = len(subs)
            for nfile in range(num_files):
                new_fname = fname[:-2] + '_subr_' + str(nfile) + '.f'
                new_fnames.append(new_fname)
                if not newer(fname, new_fname):
                    continue
                with open(new_fname, 'wb') as fn:
                    if nfile + 1 == num_files:
                        fn.writelines(lines[subs[nfile]:])
                    else:
                        fn.writelines(lines[subs[nfile]:subs[nfile+1]])

        return new_fnames

    exclude_pattern = re.compile('_subr_[0-9]')
    source_fnames = [f for f in sorted(glob.glob(os.path.join(source_dir, '*.f')))
                             if not exclude_pattern.search(os.path.basename(f))]
    fnames = []
    for source_fname in source_fnames:
        created_files = split_file(source_fname)
        if created_files is not None:
            for cfile in created_files:
                fnames.append(os.path.basename(cfile))

    return fnames
