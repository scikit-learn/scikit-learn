# -*- coding: utf-8 -*-
"""Release data for the IPython project."""

#-----------------------------------------------------------------------------
#  Copyright (c) 2008, IPython Development Team.
#  Copyright (c) 2001, Fernando Perez <fernando.perez@colorado.edu>
#  Copyright (c) 2001, Janko Hauser <jhauser@zscout.de>
#  Copyright (c) 2001, Nathaniel Gray <n8gray@caltech.edu>
#
#  Distributed under the terms of the Modified BSD License.
#
#  The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

# Name of the package for release purposes.  This is the name which labels
# the tarballs and RPMs made by distutils, so it's best to lowercase it.
name = 'ipython'

# IPython version information.  An empty _version_extra corresponds to a full
# release.  'dev' as a _version_extra string means this is a development
# version
_version_major = 6
_version_minor = 4
_version_patch = 0
_version_extra = '.dev'
# _version_extra = 'rc2'
_version_extra = ''  # Uncomment this for full releases

# Construct full version string from these.
_ver = [_version_major, _version_minor, _version_patch]

__version__ = '.'.join(map(str, _ver))
if _version_extra:
    __version__ = __version__  + _version_extra

version = __version__  # backwards compatibility name
version_info = (_version_major, _version_minor, _version_patch, _version_extra)

# Change this when incrementing the kernel protocol version
kernel_protocol_version_info = (5, 0)
kernel_protocol_version = "%i.%i" % kernel_protocol_version_info

description = "IPython: Productive Interactive Computing"

long_description = \
"""
IPython provides a rich toolkit to help you make the most out of using Python
interactively.  Its main components are:

* A powerful interactive Python shell
* A `Jupyter <http://jupyter.org/>`_ kernel to work with Python code in Jupyter
  notebooks and other interactive frontends.

The enhanced interactive Python shells have the following main features:

* Comprehensive object introspection.

* Input history, persistent across sessions.

* Caching of output results during a session with automatically generated
  references.

* Extensible tab completion, with support by default for completion of python
  variables and keywords, filenames and function keywords.

* Extensible system of 'magic' commands for controlling the environment and
  performing many tasks related either to IPython or the operating system.

* A rich configuration system with easy switching between different setups
  (simpler than changing $PYTHONSTARTUP environment variables every time).

* Session logging and reloading.

* Extensible syntax processing for special purpose situations.

* Access to the system shell with user-extensible alias system.

* Easily embeddable in other Python programs and GUIs.

* Integrated access to the pdb debugger and the Python profiler.

The latest development version is always available from IPython's `GitHub
site <http://github.com/ipython>`_.
"""

license = 'BSD'

authors = {'Fernando' : ('Fernando Perez','fperez.net@gmail.com'),
           'Janko'    : ('Janko Hauser','jhauser@zscout.de'),
           'Nathan'   : ('Nathaniel Gray','n8gray@caltech.edu'),
           'Ville'    : ('Ville Vainio','vivainio@gmail.com'),
           'Brian'    : ('Brian E Granger', 'ellisonbg@gmail.com'),
           'Min'      : ('Min Ragan-Kelley', 'benjaminrk@gmail.com'),
           'Thomas'   : ('Thomas A. Kluyver', 'takowl@gmail.com'),
           'Jorgen'   : ('Jorgen Stenarson', 'jorgen.stenarson@bostream.nu'),
           'Matthias' : ('Matthias Bussonnier', 'bussonniermatthias@gmail.com'),
           }

author = 'The IPython Development Team'

author_email = 'ipython-dev@python.org'

url = 'https://ipython.org'


platforms = ['Linux','Mac OSX','Windows']

keywords = ['Interactive','Interpreter','Shell', 'Embedding']

classifiers = [
    'Framework :: IPython',
    'Intended Audience :: Developers',
    'Intended Audience :: Science/Research',
    'License :: OSI Approved :: BSD License',
    'Programming Language :: Python',
    'Programming Language :: Python :: 3',
    'Programming Language :: Python :: 3 :: Only',
    'Topic :: System :: Shells'
    ]
