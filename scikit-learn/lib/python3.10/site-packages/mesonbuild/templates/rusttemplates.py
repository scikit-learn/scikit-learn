# SPDX-License-Identifier: Apache-2.0
# Copyright 2019 The Meson development team
# Copyright © 2023-2025 Intel Corporation

from __future__ import annotations

import typing as T

from mesonbuild.templates.sampleimpl import FileImpl

if T.TYPE_CHECKING:
    from ..minit import Arguments


lib_rust_template = '''#![crate_name = "{crate_file}"]

/* This function will not be exported and is not
 * directly callable by users of this library.
 */
fn internal_function() -> i32 {{
    return 0;
}}

pub fn {function_name}() -> i32 {{
    return internal_function();
}}

#[cfg(test)]
mod tests {{
    use super::*;

    #[test]
    fn test_function() {{
        assert_eq!({function_name}(), 0);
    }}
}}
'''


lib_rust_meson_template = '''project(
  '{project_name}',
  'rust',
  version : '{version}',
  meson_version : '>= {meson_version}',
  default_options : ['rust_std=2021', 'warning_level=3'],
)

rust = import('rust')

dependencies = [{dependencies}
]

lib = static_library(
  '{lib_name}',
  '{source_file}',
  dependencies : dependencies,
  install : true,
)

rust.test('{test_name}', lib)

# Make this library usable as a Meson subproject.
{ltoken}_dep = declare_dependency(
  include_directories : include_directories('.'),
  dependencies : dependencies,
  link_with : lib,
)
meson.override_dependency('{project_name}', {ltoken}_dep)
'''

hello_rust_template = '''
fn main() {{
    let project_name = "{project_name}";
    println!("This is project {{}}.\\n", project_name);
}}
'''

hello_rust_meson_template = '''project(
  '{project_name}',
  'rust',
  version : '{version}',
  meson_version : '>= {meson_version}',
  default_options : ['rust_std=2021', 'warning_level=3'],
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


class RustProject(FileImpl):

    source_ext = 'rs'
    exe_template = hello_rust_template
    exe_meson_template = hello_rust_meson_template
    lib_template = lib_rust_template
    lib_test_template = None
    lib_meson_template = lib_rust_meson_template

    def __init__(self, args: Arguments):
        super().__init__(args)
        self.meson_version = '1.3.0'

    def lib_kwargs(self) -> T.Dict[str, str]:
        kwargs = super().lib_kwargs()
        kwargs['crate_file'] = self.lowercase_token
        return kwargs
