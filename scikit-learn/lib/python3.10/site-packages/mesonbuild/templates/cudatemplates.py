# SPDX-License-Identifier: Apache-2.0
# Copyright 2019 The Meson development team
# Copyright © 2023-2025 Intel Corporation

from __future__ import annotations
import typing as T

from mesonbuild.templates.sampleimpl import FileHeaderImpl

if T.TYPE_CHECKING:
    from ..minit import Arguments


hello_cuda_template = '''#include <iostream>

#define PROJECT_NAME "{project_name}"

int main(int argc, char **argv) {{
    if (argc != 1) {{
        std::cout << argv[0] << " takes no arguments.\\n";
        return 1;
    }}
    std::cout << "This is project " << PROJECT_NAME << ".\\n";
    return 0;
}}
'''

hello_cuda_meson_template = '''project(
  '{project_name}',
  ['cuda', 'cpp'],
  version : '{version}',
  meson_version : '>= {meson_version}',
  default_options : ['warning_level=3', 'cpp_std=c++14'],
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

lib_h_template = '''#pragma once
#if defined _WIN32 || defined __CYGWIN__
  #ifdef BUILDING_{utoken}
    #define {utoken}_PUBLIC __declspec(dllexport)
  #else
    #define {utoken}_PUBLIC __declspec(dllimport)
  #endif
#else
  #ifdef BUILDING_{utoken}
      #define {utoken}_PUBLIC __attribute__ ((visibility ("default")))
  #else
      #define {utoken}_PUBLIC
  #endif
#endif

namespace {namespace} {{

class {utoken}_PUBLIC {class_name} {{

public:
  {class_name}();
  int get_number() const;

private:

  int number;

}};

}}

'''

lib_cuda_template = '''#include <{header_file}>

namespace {namespace} {{

{class_name}::{class_name}() {{
    number = 6;
}}

int {class_name}::get_number() const {{
  return number;
}}

}}
'''

lib_cuda_test_template = '''#include <{header_file}>
#include <iostream>

int main(int argc, char **argv) {{
    if (argc != 1) {{
        std::cout << argv[0] << " takes no arguments.\\n";
        return 1;
    }}
    {namespace}::{class_name} c;
    return c.get_number() != 6;
}}
'''

lib_cuda_meson_template = '''project(
  '{project_name}',
  ['cuda', 'cpp'],
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
  cpp_shared_args : lib_args,
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

# Make this library usable from the system's
# package manager.
install_headers('{header_file}', subdir : '{header_dir}')

pkg_mod = import('pkgconfig')
pkg_mod.generate(
  lib,
  description : 'Meson sample project.',
  subdirs : '{header_dir}',
)
'''


class CudaProject(FileHeaderImpl):

    source_ext = 'cu'
    header_ext = 'h'
    exe_template = hello_cuda_template
    exe_meson_template = hello_cuda_meson_template
    lib_template = lib_cuda_template
    lib_header_template = lib_h_template
    lib_test_template = lib_cuda_test_template
    lib_meson_template = lib_cuda_meson_template

    def __init__(self, args: Arguments):
        super().__init__(args)
        self.meson_version = '1.3.0'
