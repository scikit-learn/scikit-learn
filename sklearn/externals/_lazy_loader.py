# BSD 3-Clause License

# Copyright (c) 2022--2023, Scientific Python project
# All rights reserved.

# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:

# 1. Redistributions of source code must retain the above copyright notice, this
#    list of conditions and the following disclaimer.

# 2. Redistributions in binary form must reproduce the above copyright notice,
#    this list of conditions and the following disclaimer in the documentation
#    and/or other materials provided with the distribution.

# 3. Neither the name of the copyright holder nor the names of its
#    contributors may be used to endorse or promote products derived from
#    this software without specific prior written permission.

# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
# FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
# DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
# SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
# CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
# OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

# lazy_loader
# ===========

# Makes it easy to load subpackages and functions on demand.

import importlib
import importlib.util
import os
import sys

__all__ = ["attach"]


def attach(package_name, submodules=None, submod_attrs=None):
    """Attach lazily loaded submodules, functions, or other attributes.

    Typically, modules import submodules and attributes as follows::

      import mysubmodule
      import anothersubmodule

      from .foo import someattr

    The idea is to replace a package's `__getattr__`, `__dir__`, and
    `__all__`, such that all imports work exactly the way they would
    with normal imports, except that the import occurs upon first use.

    The typical way to call this function, replacing the above imports, is::

      __getattr__, __dir__, __all__ = lazy.attach(
        __name__,
        ['mysubmodule', 'anothersubmodule'],
        {'foo': ['someattr']}
      )

    This functionality requires Python 3.7 or higher.

    Parameters
    ----------
    package_name : str
        Typically use ``__name__``.
    submodules : set
        List of submodules to attach.
    submod_attrs : dict
        Dictionary of submodule -> list of attributes / functions.
        These attributes are imported as they are used.

    Returns
    -------
    __getattr__, __dir__, __all__

    """
    if submod_attrs is None:
        submod_attrs = {}

    if submodules is None:
        submodules = set()
    else:
        submodules = set(submodules)

    attr_to_modules = {
        attr: mod for mod, attrs in submod_attrs.items() for attr in attrs
    }

    __all__ = list(sorted(submodules | attr_to_modules.keys()))

    def __getattr__(name):
        if name in submodules:
            return importlib.import_module(f"{package_name}.{name}")
        elif name in attr_to_modules:
            submod_path = f"{package_name}.{attr_to_modules[name]}"
            submod = importlib.import_module(submod_path)
            attr = getattr(submod, name)

            # If the attribute lives in a file (module) with the same
            # name as the attribute, ensure that the attribute and *not*
            # the module is accessible on the package.
            if name == attr_to_modules[name]:
                pkg = sys.modules[package_name]
                pkg.__dict__[name] = attr

            return attr
        else:
            raise AttributeError(f"No {package_name} attribute {name}")

    def __dir__():
        return __all__

    if os.environ.get("EAGER_IMPORT", ""):
        for attr in set(attr_to_modules.keys()) | submodules:
            __getattr__(attr)

    return __getattr__, __dir__, list(__all__)
