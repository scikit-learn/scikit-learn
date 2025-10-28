# SPDX-License-Identifier: Apache-2.0
# Copyright 2019 The Meson development team
# Copyright © 2023-2025 Intel Corporation

from __future__ import annotations

from mesonbuild.templates.sampleimpl import ClassImpl


hello_java_template = '''

public class {class_name} {{
    final static String PROJECT_NAME = "{project_name}";

    public static void main (String args[]) {{
        if (args.length != 0) {{
            System.out.println(args + " takes no arguments.");
            System.exit(0);
        }}
        System.out.println("This is project " + PROJECT_NAME + ".");
        System.exit(0);
    }}
}}

'''

hello_java_meson_template = '''project(
  '{project_name}',
  'java',
  version : '{version}',
  meson_version : '>= {meson_version}',
  default_options : ['warning_level=3'],
)

dependencies = [{dependencies}
]

exe = jar(
  '{exe_name}',
  '{source_name}',
  main_class : '{exe_name}',
  dependencies : dependencies,
  install : true,
)

test('basic', exe)
'''

lib_java_template = '''

public class {class_name} {{
    final static int number = 6;

    public final int get_number() {{
      return number;
    }}
}}

'''

lib_java_test_template = '''

public class {class_test} {{
    public static void main (String args[]) {{
        if (args.length != 0) {{
            System.out.println(args + " takes no arguments.");
            System.exit(1);
        }}

        {class_name} c = new {class_name}();
        Boolean result = true;
        System.exit(result.compareTo(c.get_number() != 6));
    }}
}}

'''

lib_java_meson_template = '''project(
  '{project_name}',
  'java',
  version : '{version}',
  meson_version : '>= {meson_version}',
  default_options : ['warning_level=3'],
)

dependencies = [{dependencies}
]

jarlib = jar(
  '{class_name}',
  '{source_file}',
  dependencies : dependencies,
  main_class : '{class_name}',
  install : true,
)

test_jar = jar(
  '{class_test}',
  '{test_source_file}',
  main_class : '{class_test}',
  dependencies : dependencies,
  link_with : jarlib,
)
test('{test_name}', test_jar)

# Make this library usable as a Meson subproject.
{ltoken}_dep = declare_dependency(
  include_directories : include_directories('.'),
  dependencies : dependencies,
  link_with : jarlib,
)
meson.override_dependency('{project_name}', {ltoken}_dep)
'''


class JavaProject(ClassImpl):

    source_ext = 'java'
    exe_template = hello_java_template
    exe_meson_template = hello_java_meson_template
    lib_template = lib_java_template
    lib_test_template = lib_java_test_template
    lib_meson_template = lib_java_meson_template
