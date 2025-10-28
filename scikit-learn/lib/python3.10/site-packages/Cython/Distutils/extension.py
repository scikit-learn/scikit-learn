"""Pyrex.Distutils.extension

Provides a modified Extension class, that understands how to describe
Pyrex extension modules in setup scripts."""

__revision__ = "$Id:$"

import distutils.extension as _Extension


class Extension(_Extension.Extension):
    # When adding arguments to this constructor, be sure to update
    # user_options.extend in build_ext.py.
    def __init__(self, name, sources,
                 include_dirs=None,
                 define_macros=None,
                 undef_macros=None,
                 library_dirs=None,
                 libraries=None,
                 runtime_library_dirs=None,
                 extra_objects=None,
                 extra_compile_args=None,
                 extra_link_args=None,
                 export_symbols=None,
                 #swig_opts=None,
                 depends=None,
                 language=None,
                 cython_include_dirs=None,
                 cython_directives=None,
                 cython_create_listing=False,
                 cython_line_directives=False,
                 cython_cplus=False,
                 cython_c_in_temp=False,
                 cython_gen_pxi=False,
                 cython_gdb=False,
                 no_c_in_traceback=False,
                 cython_compile_time_env=None,
                 **kw):

        # Translate pyrex_X to cython_X for backwards compatibility.
        had_pyrex_options = False
        for key in list(kw):
            if key.startswith('pyrex_'):
                had_pyrex_options = True
                kw['cython' + key[5:]] = kw.pop(key)
        if had_pyrex_options:
            Extension.__init__(
                self, name, sources,
                include_dirs=include_dirs,
                define_macros=define_macros,
                undef_macros=undef_macros,
                library_dirs=library_dirs,
                libraries=libraries,
                runtime_library_dirs=runtime_library_dirs,
                extra_objects=extra_objects,
                extra_compile_args=extra_compile_args,
                extra_link_args=extra_link_args,
                export_symbols=export_symbols,
                #swig_opts=swig_opts,
                depends=depends,
                language=language,
                no_c_in_traceback=no_c_in_traceback,
                **kw)
            return

        _Extension.Extension.__init__(
            self, name, sources,
            include_dirs=include_dirs,
            define_macros=define_macros,
            undef_macros=undef_macros,
            library_dirs=library_dirs,
            libraries=libraries,
            runtime_library_dirs=runtime_library_dirs,
            extra_objects=extra_objects,
            extra_compile_args=extra_compile_args,
            extra_link_args=extra_link_args,
            export_symbols=export_symbols,
            #swig_opts=swig_opts,
            depends=depends,
            language=language,
            **kw)

        self.cython_include_dirs = cython_include_dirs or []
        self.cython_directives = cython_directives or {}
        self.cython_create_listing = cython_create_listing
        self.cython_line_directives = cython_line_directives
        self.cython_cplus = cython_cplus
        self.cython_c_in_temp = cython_c_in_temp
        self.cython_gen_pxi = cython_gen_pxi
        self.cython_gdb = cython_gdb
        self.no_c_in_traceback = no_c_in_traceback
        self.cython_compile_time_env = cython_compile_time_env

# class Extension

read_setup_file = _Extension.read_setup_file
