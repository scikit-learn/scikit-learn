"""
Pythran integration into IPython.

* provides the %%pythran magic function to ipython
"""
# -----------------------------------------------------------------------------
# Copyright (C) 2010-2011, IPython Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# -----------------------------------------------------------------------------

import hashlib
import importlib
from IPython.core.magic import Magics, magics_class, cell_magic
from IPython.core import magic_arguments

import pythran


@magics_class
class PythranMagics(Magics):

    """ Class to make it possible to use pythran as a magic IPython command."""

    def __init__(self, shell):
        """ Init the pythran magic stuff. """
        super(PythranMagics, self).__init__(shell)
        self._reloads = {}

    def _import_all(self, module):
        """ Import only globals modules. """
        self.shell.push({k: v for k, v in module.__dict__.items()
                         if not k.startswith('__')})

    @magic_arguments.magic_arguments()
    @magic_arguments.argument('-D', action='append', default=[])
    @magic_arguments.argument('-O', action='append', default=[])
    @magic_arguments.argument('-m', action='append', default=[])
    @magic_arguments.argument('-W', action='append', default=[])
    @magic_arguments.argument('-f', action='append', default=[])
    @cell_magic
    def pythran(self, line, cell):
        """
        Compile and import everything from a Pythran code cell.

        %%pythran
        #pythran export foo(int)
        def foo(x):
            return x + x
        """
        args = magic_arguments.parse_argstring(self.pythran, line)
        kwargs = {}
        if args.D:
            kwargs['define_macros'] = args.D
        for v in "OmWf":
            args_v = getattr(args, v)
            for target in ('extra_compile_args', 'extra_link_args'):
                kwargs.setdefault(target, []).extend('-{}{}'.format(v, x)
                                                     for x in args_v)

        m = hashlib.md5()
        m.update(line.encode('utf-8'))
        m.update(cell.encode('utf-8'))
        module_name = "pythranized_" + m.hexdigest()
        module_path = pythran.compile_pythrancode(module_name, cell, **kwargs)
        loader = importlib.machinery.ExtensionFileLoader(module_name, module_path)
        spec = importlib.machinery.ModuleSpec(name=module_name, loader=loader,
                                              origin=module_path)
        module = importlib._bootstrap._load(spec)
        self._import_all(module)


def load_ipython_extension(ipython):
    """Load the extension in IPython."""
    ipython.register_magics(PythranMagics)
