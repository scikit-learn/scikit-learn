"""Utility functions for finding modules

Utility functions for finding modules on sys.path.

`find_mod` finds named module on sys.path.

`get_init` helper function that finds __init__ file in a directory.

`find_module` variant of imp.find_module in std_lib that only returns
path to module and not an open file object as well.



"""
#-----------------------------------------------------------------------------
# Copyright (c) 2011, the IPython Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
# Imports
#-----------------------------------------------------------------------------

# Stdlib imports
import imp
import os

# Third-party imports

# Our own imports


#-----------------------------------------------------------------------------
# Globals and constants
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
# Local utilities
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
# Classes and functions
#-----------------------------------------------------------------------------
def find_module(name, path=None):
    """imp.find_module variant that only return path of module.
    
    The `imp.find_module` returns a filehandle that we are not interested in.
    Also we ignore any bytecode files that `imp.find_module` finds.

    Parameters
    ----------
    name : str
        name of module to locate
    path : list of str
        list of paths to search for `name`. If path=None then search sys.path

    Returns
    -------
    filename : str
        Return full path of module or None if module is missing or does not have
        .py or .pyw extension
    """
    if name is None:
        return None
    try:
        file, filename, _ = imp.find_module(name, path)
    except ImportError:
        return None
    if file is None:
        return filename
    else:
        file.close()
    if os.path.splitext(filename)[1] in [".py", ".pyc"]:
        return filename
    else:
        return None

def get_init(dirname):
    """Get __init__ file path for module directory
    
    Parameters
    ----------
    dirname : str
        Find the __init__ file in directory `dirname`

    Returns
    -------
    init_path : str
        Path to __init__ file
    """
    fbase =  os.path.join(dirname, "__init__")
    for ext in [".py", ".pyw"]:
        fname = fbase + ext
        if os.path.isfile(fname):
            return fname


def find_mod(module_name):
    """Find module `module_name` on sys.path
    
    Return the path to module `module_name`. If `module_name` refers to
    a module directory then return path to __init__ file. Return full 
    path of module or None if module is missing or does not have .py or .pyw
    extension. We are not interested in running bytecode.
    
    Parameters
    ----------
    module_name : str
    
    Returns
    -------
    modulepath : str
        Path to module `module_name`.
    """
    parts = module_name.split(".")
    basepath = find_module(parts[0])
    for submodname in parts[1:]:
        basepath = find_module(submodname, [basepath])
    if basepath and os.path.isdir(basepath):
        basepath = get_init(basepath)
    return basepath
