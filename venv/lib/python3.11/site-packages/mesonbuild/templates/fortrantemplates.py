# SPDX-License-Identifier: Apache-2.0
# Copyright 2019 The Meson development team

from __future__ import annotations

from mesonbuild.templates.sampleimpl import FileImpl

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

lib_fortran_meson_template = '''project('{project_name}', 'fortran',
  version : '{version}',
  default_options : ['warning_level=3'])

# These arguments are only used to build the shared library
# not the executables that use the library.
lib_args = ['-DBUILDING_{utoken}']

shlib = shared_library('{lib_name}', '{source_file}',
  install : true,
  fortran_args : lib_args,
  gnu_symbol_visibility : 'hidden',
)

test_exe = executable('{test_exe_name}', '{test_source_file}',
  link_with : shlib)
test('{test_name}', test_exe)

# Make this library usable as a Meson subproject.
{ltoken}_dep = declare_dependency(
  include_directories: include_directories('.'),
  link_with : shlib)

pkg_mod = import('pkgconfig')
pkg_mod.generate(
  name : '{project_name}',
  filebase : '{ltoken}',
  description : 'Meson sample project.',
  subdirs : '{header_dir}',
  libraries : shlib,
  version : '{version}',
)
'''

hello_fortran_template = '''
implicit none

character(len=*), parameter :: PROJECT_NAME = "{project_name}"

print *,"This is project ", PROJECT_NAME

end program
'''

hello_fortran_meson_template = '''project('{project_name}', 'fortran',
  version : '{version}',
  default_options : ['warning_level=3'])

exe = executable('{exe_name}', '{source_name}',
  install : true)

test('basic', exe)
'''


class FortranProject(FileImpl):

    source_ext = 'f90'
    exe_template = hello_fortran_template
    exe_meson_template = hello_fortran_meson_template
    lib_template = lib_fortran_template
    lib_meson_template = lib_fortran_meson_template
    lib_test_template = lib_fortran_test_template
