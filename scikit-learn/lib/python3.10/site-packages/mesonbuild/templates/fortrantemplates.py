# SPDX-License-Identifier: Apache-2.0
# Copyright 2019 The Meson development team
# Copyright © 2023-2025 Intel Corporation

from __future__ import annotations
import typing as T

from mesonbuild.templates.sampleimpl import FileImpl

if T.TYPE_CHECKING:
    from ..minit import Arguments

lib_fortran_template = '''
! This procedure will not be exported and is not
! directly callable by users of this library.

module modfoo

implicit none
private
public :: {function_name}

contains

integer function internal_function()
    internal_function = 0
end function internal_function

integer function {function_name}()
    {function_name} = internal_function()
end function {function_name}

end module modfoo
'''

lib_fortran_test_template = '''
use modfoo

print *,{function_name}()

end program
'''

lib_fortran_meson_template = '''project(
  '{project_name}',
  'fortran',
  version : '{version}',
  meson_version : '>= {meson_version}',
  default_options : ['warning_level=3'],
)

# These arguments are only used to build the shared library
# not the executables that use the library.
lib_args = ['-DBUILDING_{utoken}']

dependencies = [{dependencies}
]

lib = library(
  '{lib_name}',
  '{source_file}',
  install : true,
  fortran_shared_args : lib_args,
  gnu_symbol_visibility : 'hidden',
  dependencies : dependencies,
)

test_exe = executable(
  '{test_exe_name}',
  '{test_source_file}',
  link_with : lib,
  dependencies : dependencies,
)
test('{test_name}', test_exe)

# Make this library usable as a Meson subproject.
{ltoken}_dep = declare_dependency(
  include_directories : include_directories('.'),
  dependencies : dependencies,
  link_with : lib,
)
meson.override_dependency('{project_name}', {ltoken}_dep)

pkg_mod = import('pkgconfig')
pkg_mod.generate(
  lib,
  description : 'Meson sample project.',
  subdirs : '{header_dir}',
)
'''

hello_fortran_template = '''
implicit none

character(len=*), parameter :: PROJECT_NAME = "{project_name}"

print *,"This is project ", PROJECT_NAME

end program
'''

hello_fortran_meson_template = '''project(
  '{project_name}',
  'fortran',
  version : '{version}',
  meson_version : '>= {meson_version}',
  default_options : ['warning_level=3'],
)

dependencies = [{dependencies}
]

exe = executable(
  '{exe_name}',
  '{source_name}',
  dependencies : dependencies,
  install : true,
)

test('basic', exe)
'''


class FortranProject(FileImpl):

    source_ext = 'f90'
    exe_template = hello_fortran_template
    exe_meson_template = hello_fortran_meson_template
    lib_template = lib_fortran_template
    lib_meson_template = lib_fortran_meson_template
    lib_test_template = lib_fortran_test_template

    def __init__(self, args: Arguments):
        super().__init__(args)
        self.meson_version = '1.3.0'
