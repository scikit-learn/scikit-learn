# These functions are based on the following code:
# https://git.savannah.gnu.org/gitweb/?p=autoconf-archive.git;a=blob_plain;f=m4/ax_gcc_func_attribute.m4,
# which is licensed under the following terms:
#
#   Copyright (c) 2013 Gabriele Svelto <gabriele.svelto@gmail.com>
#
#   Copying and distribution of this file, with or without modification, are
#   permitted in any medium without royalty provided the copyright notice
#   and this notice are preserved.  This file is offered as-is, without any
#   warranty.
#

C_FUNC_ATTRIBUTES = {
    'alias': '''
        int foo(void) { return 0; }
        int bar(void) __attribute__((alias("foo")));''',
    'aligned':
        'int foo(void) __attribute__((aligned(32)));',
    'alloc_size':
        'void *foo(int a) __attribute__((alloc_size(1)));',
    'always_inline':
        'inline __attribute__((always_inline)) int foo(void) { return 0; }',
    'artificial':
        'inline __attribute__((artificial)) int foo(void) { return 0; }',
    'cold':
        'int foo(void) __attribute__((cold));',
    'const':
        'int foo(void) __attribute__((const));',
    'constructor':
        'int foo(void) __attribute__((constructor));',
    'constructor_priority':
        'int foo( void ) __attribute__((__constructor__(65535/2)));',
    'deprecated':
        'int foo(void) __attribute__((deprecated("")));',
    'destructor':
        'int foo(void) __attribute__((destructor));',
    'dllexport':
        '__declspec(dllexport) int foo(void) { return 0; }',
    'dllimport':
        '__declspec(dllimport) int foo(void);',
    'error':
        'int foo(void) __attribute__((error("")));',
    'externally_visible':
        'int foo(void) __attribute__((externally_visible));',
    'fallthrough': '''
        int foo( void ) {
          switch (0) {
            case 1: __attribute__((fallthrough));
            case 2: break;
          }
          return 0;
        };''',
    'flatten':
        'int foo(void) __attribute__((flatten));',
    'format':
        'int foo(const char * p, ...) __attribute__((format(printf, 1, 2)));',
    'format_arg':
        'char * foo(const char * p) __attribute__((format_arg(1)));',
    'force_align_arg_pointer':
        '__attribute__((force_align_arg_pointer)) int foo(void) { return 0; }',
    'gnu_inline':
        'inline __attribute__((gnu_inline)) int foo(void) { return 0; }',
    'hot':
        'int foo(void) __attribute__((hot));',
    'ifunc':
        ('int my_foo(void) { return 0; }'
         'static int (*resolve_foo(void))(void) { return my_foo; }'
         'int foo(void) __attribute__((ifunc("resolve_foo")));'),
    'leaf':
        '__attribute__((leaf)) int foo(void) { return 0; }',
    'malloc':
        'int *foo(void) __attribute__((malloc));',
    'noclone':
        'int foo(void) __attribute__((noclone));',
    'noinline':
        '__attribute__((noinline)) int foo(void) { return 0; }',
    'nonnull':
        'int foo(char * p) __attribute__((nonnull(1)));',
    'noreturn':
        'int foo(void) __attribute__((noreturn));',
    'nothrow':
        'int foo(void) __attribute__((nothrow));',
    'null_terminated_string_arg':
        'int foo(const char * p) __attribute__((null_terminated_string_arg(1)));',
    'optimize':
        '__attribute__((optimize(3))) int foo(void) { return 0; }',
    'packed':
        'struct __attribute__((packed)) foo { int bar; };',
    'pure':
        'int foo(void) __attribute__((pure));',
    'returns_nonnull':
        'int *foo(void) __attribute__((returns_nonnull));',
    'section': '''
        #if defined(__APPLE__) && defined(__MACH__)
            extern int foo __attribute__((section("__BAR,__bar")));
        #else
            extern int foo __attribute__((section(".bar")));
        #endif''',
    'sentinel':
        'int foo(const char *bar, ...) __attribute__((sentinel));',
    'unused':
        'int foo(void) __attribute__((unused));',
    'used':
        'int foo(void) __attribute__((used));',
    'vector_size':
        '__attribute__((vector_size(32))); int foo(void) { return 0; }',
    'visibility': '''
        int foo_def(void) __attribute__((visibility("default"))); int foo_def(void) { return 0; }
        int foo_hid(void) __attribute__((visibility("hidden"))); int foo_hid(void) { return 0; }
        int foo_int(void) __attribute__((visibility("internal"))); int foo_int(void) { return 0; }''',
    'visibility:default':
        'int foo(void) __attribute__((visibility("default"))); int foo(void) { return 0; }',
    'visibility:hidden':
        'int foo(void) __attribute__((visibility("hidden"))); int foo(void) { return 0; }',
    'visibility:internal':
        'int foo(void) __attribute__((visibility("internal"))); int foo(void) { return 0; }',
    'visibility:protected':
        'int foo(void) __attribute__((visibility("protected"))); int foo(void) { return 0; }',
    'warning':
        'int foo(void) __attribute__((warning("")));',
    'warn_unused_result':
        'int foo(void) __attribute__((warn_unused_result));',
    'weak':
        'int foo(void) __attribute__((weak));',
    'weakref': '''
        static int foo(void) { return 0; }
        static int var(void) __attribute__((weakref("foo")));''',
    'retain': '__attribute__((retain)) int x;',
}

CXX_FUNC_ATTRIBUTES = {
    # Alias must be applied to the mangled name in C++
    'alias':
        ('extern "C" {'
         'int foo(void) { return 0; }'
         '}'
         'int bar(void) __attribute__((alias("foo")));'
         ),
    'ifunc':
        ('extern "C" {'
         'int my_foo(void) { return 0; }'
         'static int (*resolve_foo(void))(void) { return my_foo; }'
         '}'
         'int foo(void) __attribute__((ifunc("resolve_foo")));'),
}
