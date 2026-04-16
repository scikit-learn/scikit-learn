# SPDX-License-Identifier: Apache-2.0
# Copyright 2019 The Meson development team
# Copyright © 2023-2025 Intel Corporation

from __future__ import annotations

from mesonbuild.templates.sampleimpl import ClassImpl


hello_cs_template = '''using System;

public class {class_name} {{
    const String PROJECT_NAME = "{project_name}";

    static int Main(String[] args) {{
      if (args.Length > 0) {{
          System.Console.WriteLine(String.Format("{project_name} takes no arguments.."));
          return 1;
      }}
      Console.WriteLine(String.Format("This is project {{0}}.", PROJECT_NAME));
      return 0;
    }}
}}

'''

hello_cs_meson_template = '''project(
  '{project_name}',
  'cs',
  version : '{version}',
  meson_version : '>= {meson_version}',
  default_options : ['warning_level=3'],
)

dependencies = [{dependencies}
]

sources = [{source_files}
]

exe = executable(
  '{exe_name}',
  sources,
  install : true,
  dependencies : dependencies,
)

test('basic', exe)
'''

lib_cs_template = '''
public class {class_name} {{
    private const int number = 6;

    public int get_number() {{
      return number;
    }}
}}

'''

lib_cs_test_template = '''using System;

public class {class_test} {{
    static int Main(String[] args) {{
      if (args.Length > 0) {{
          System.Console.WriteLine("{project_name} takes no arguments..");
          return 1;
      }}
      {class_name} c = new {class_name}();
      Boolean result = true;
      return result.CompareTo(c.get_number() != 6);
    }}
}}

'''

lib_cs_meson_template = '''project(
  '{project_name}',
  'cs',
  version : '{version}',
  meson_version : '>= {meson_version}',
  default_options : ['warning_level=3'],
)

dependencies = [{dependencies}
]

sources = [{source_files}
]

stlib = shared_library(
  '{lib_name}',
  sources,
  dependencies : dependencies,
  install : true,
)

test_exe = executable(
  '{test_exe_name}',
  '{test_source_file}',
  dependencies : dependencies,
  link_with : stlib,
)
test('{test_name}', test_exe)

# Make this library usable as a Meson subproject.
{ltoken}_dep = declare_dependency(
  include_directories : include_directories('.'),
  dependencies : dependencies,
  link_with : stlib,
)
meson.override_dependency('{project_name}', {ltoken}_dep)

'''


class CSharpProject(ClassImpl):

    source_ext = 'cs'
    exe_template = hello_cs_template
    exe_meson_template = hello_cs_meson_template
    lib_template = lib_cs_template
    lib_test_template = lib_cs_test_template
    lib_meson_template = lib_cs_meson_template
