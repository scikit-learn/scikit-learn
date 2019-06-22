# coding: utf-8

from __future__ import print_function, absolute_import, division, unicode_literals

if False:  # MYPY
    from typing import Dict, Any  # NOQA

_package_data = dict(
    full_package_name='ruamel.yaml',
    version_info=(0, 15, 97),
    __version__='0.15.97',
    author='Anthon van der Neut',
    author_email='a.van.der.neut@ruamel.eu',
    description='ruamel.yaml is a YAML parser/emitter that supports roundtrip preservation of comments, seq/map flow style, and map key order',  # NOQA
    entry_points=None,
    since=2014,
    extras_require={':platform_python_implementation=="CPython" and python_version<="2.7"': [
            'ruamel.ordereddict',
        ], 'jinja2': ['ruamel.yaml.jinja2>=0.2'], 'docs': ['ryd']},
    ext_modules=[
            dict(
                name='_ruamel_yaml',
                src=[
                        'ext/_ruamel_yaml.c',
                        'ext/api.c',
                        'ext/writer.c',
                        'ext/dumper.c',
                        'ext/loader.c',
                        'ext/reader.c',
                        'ext/scanner.c',
                        'ext/parser.c',
                        'ext/emitter.c',
                ],
                lib=[],
                test="""
            int main(int argc, char* argv[])
            {
              /* prevent warning */
              return 0;
            }
            """,
            ),
    ],
    # NOQA
    # test='#include "ext/yaml.h"\n\nint main(int argc, char* argv[])\n{\nyaml_parser_t parser;\nparser = parser;  /* prevent warning */\nreturn 0;\n}\n',  # NOQA
    classifiers=[
            'Programming Language :: Python :: 2.7',
            'Programming Language :: Python :: 3.5',
            'Programming Language :: Python :: 3.6',
            'Programming Language :: Python :: 3.7',
            'Programming Language :: Python :: Implementation :: CPython',
            'Programming Language :: Python :: Implementation :: PyPy',
            'Programming Language :: Python :: Implementation :: Jython',
            'Topic :: Software Development :: Libraries :: Python Modules',
            'Topic :: Text Processing :: Markup',
    ],
    keywords='yaml 1.2 parser round-trip preserve quotes order config',
    wheels=dict(
        windows='appveyor',
        linux='libyaml-devel',
        macos='builder@macos',
    ),
    read_the_docs='yaml',
    supported=[(2, 7), (3, 5)],  # minimum
    tox=dict(
        env='*pn',  # also test narrow Python 2.7.15 for unicode patterns
        deps='ruamel.std.pathlib',
        fl8excl='_test/lib',
    ),
    rtfd='yaml',
)  # type: Dict[Any, Any]


version_info = _package_data['version_info']
__version__ = _package_data['__version__']

try:
    from .cyaml import *  # NOQA

    __with_libyaml__ = True
except (ImportError, ValueError):  # for Jython
    __with_libyaml__ = False

from ruamel.yaml.main import *  # NOQA
