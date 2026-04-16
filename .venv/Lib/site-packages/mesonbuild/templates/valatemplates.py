# SPDX-License-Identifier: Apache-2.0
# Copyright 2019 The Meson development team
# Copyright © 2023-2025 Intel Corporation

from __future__ import annotations

from mesonbuild.templates.sampleimpl import FileImpl


hello_vala_template = '''void main (string[] args) {{
    stdout.printf ("Hello {project_name}!\\n");
}}
'''

hello_vala_meson_template = '''project(
  '{project_name}',
  'vala',
  meson_version : '>= {meson_version}',
  version : '{version}',
)

dependencies = [
  dependency('glib-2.0'),
  dependency('gobject-2.0'),{dependencies}
]

sources = [{source_files}
]

exe = executable(
  '{exe_name}',
  sources,
  dependencies : dependencies,
  install : true,
)

test('basic', exe)
'''


lib_vala_template = '''namespace {namespace} {{
    public int sum(int a, int b) {{
        return(a + b);
    }}

    public int square(int a) {{
        return(a * a);
    }}
}}
'''

lib_vala_test_template = '''using {namespace};

public void main() {{
    stdout.printf("\nTesting shlib");
    stdout.printf("\n\t2 + 3 is %d", sum(2, 3));
    stdout.printf("\n\t8 squared is %d\\n", square(8));
}}
'''

lib_vala_meson_template = '''project(
  '{project_name}',
  'vala',
  meson_version : '>= {meson_version}',
  version : '{version}',
)

dependencies = [
  dependency('glib-2.0'),
  dependency('gobject-2.0'),{dependencies}
]

sources = [{source_files}
]

# These arguments are only used to build the shared library
# not the executables that use the library.
lib = shared_library(
  'foo',
  sources,
  dependencies : dependencies,
  install : true,
  install_dir : [true, true, true],
)

test_exe = executable(
  '{test_exe_name}',
  '{test_source_file}',
  dependencies : dependencies,
  link_with : lib,
)
test('{test_name}', test_exe)

# Make this library usable as a Meson subproject.
{ltoken}_dep = declare_dependency(
  include_directories : include_directories('.'),
  dependencies : dependencies,
  link_with : lib,
)
meson.override_dependency('{project_name}', {ltoken}_dep)

'''


class ValaProject(FileImpl):

    source_ext = 'vala'
    exe_template = hello_vala_template
    exe_meson_template = hello_vala_meson_template
    lib_template = lib_vala_template
    lib_test_template = lib_vala_test_template
    lib_meson_template = lib_vala_meson_template
