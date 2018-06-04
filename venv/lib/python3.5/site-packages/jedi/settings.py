"""
This module contains variables with global |jedi| settings. To change the
behavior of |jedi|, change the variables defined in :mod:`jedi.settings`.

Plugins should expose an interface so that the user can adjust the
configuration.


Example usage::

    from jedi import settings
    settings.case_insensitive_completion = True


Completion output
~~~~~~~~~~~~~~~~~

.. autodata:: case_insensitive_completion
.. autodata:: add_bracket_after_function
.. autodata:: no_completion_duplicates


Filesystem cache
~~~~~~~~~~~~~~~~

.. autodata:: cache_directory
.. autodata:: use_filesystem_cache


Parser
~~~~~~

.. autodata:: fast_parser


Dynamic stuff
~~~~~~~~~~~~~

.. autodata:: dynamic_array_additions
.. autodata:: dynamic_params
.. autodata:: dynamic_params_for_other_modules
.. autodata:: additional_dynamic_modules
.. autodata:: auto_import_modules


Caching
~~~~~~~

.. autodata:: call_signatures_validity


"""
import os
import platform

# ----------------
# completion output settings
# ----------------

case_insensitive_completion = True
"""
The completion is by default case insensitive.
"""

add_bracket_after_function = False
"""
Adds an opening bracket after a function, because that's normal behaviour.
Removed it again, because in VIM that is not very practical.
"""

no_completion_duplicates = True
"""
If set, completions with the same name don't appear in the output anymore,
but are in the `same_name_completions` attribute.
"""

# ----------------
# Filesystem cache
# ----------------

use_filesystem_cache = True
"""
Use filesystem cache to save once parsed files with pickle.
"""

if platform.system().lower() == 'windows':
    _cache_directory = os.path.join(os.getenv('APPDATA') or '~', 'Jedi',
                                    'Jedi')
elif platform.system().lower() == 'darwin':
    _cache_directory = os.path.join('~', 'Library', 'Caches', 'Jedi')
else:
    _cache_directory = os.path.join(os.getenv('XDG_CACHE_HOME') or '~/.cache',
                                    'jedi')
cache_directory = os.path.expanduser(_cache_directory)
"""
The path where the cache is stored.

On Linux, this defaults to ``~/.cache/jedi/``, on OS X to
``~/Library/Caches/Jedi/`` and on Windows to ``%APPDATA%\\Jedi\\Jedi\\``.
On Linux, if environment variable ``$XDG_CACHE_HOME`` is set,
``$XDG_CACHE_HOME/jedi`` is used instead of the default one.
"""

# ----------------
# parser
# ----------------

fast_parser = True
"""
Use the fast parser. This means that reparsing is only being done if
something has been changed e.g. to a function. If this happens, only the
function is being reparsed.
"""

# ----------------
# dynamic stuff
# ----------------

dynamic_array_additions = True
"""
check for `append`, etc. on arrays: [], {}, () as well as list/set calls.
"""

dynamic_params = True
"""
A dynamic param completion, finds the callees of the function, which define
the params of a function.
"""

dynamic_params_for_other_modules = True
"""
Do the same for other modules.
"""

additional_dynamic_modules = []
"""
Additional modules in which |jedi| checks if statements are to be found. This
is practical for IDEs, that want to administrate their modules themselves.
"""

dynamic_flow_information = True
"""
Check for `isinstance` and other information to infer a type.
"""

auto_import_modules = [
    'hashlib',  # setattr
]
"""
Modules that are not analyzed but imported, although they contain Python code.
This improves autocompletion for libraries that use ``setattr`` or
``globals()`` modifications a lot.
"""

# ----------------
# caching validity (time)
# ----------------

call_signatures_validity = 3.0
"""
Finding function calls might be slow (0.1-0.5s). This is not acceptible for
normal writing. Therefore cache it for a short time.
"""
