import tempfile
import os
import shutil

from Cython.Compiler import (
    MemoryView, Code, Options, Pipeline, Errors, Main, Symtab
)
from Cython.Compiler.StringEncoding import EncodedString
from Cython.Compiler.Scanning import FileSourceDescriptor


def create_shared_library_pipeline(context, scope, options, result):

    parse = Pipeline.parse_stage_factory(context)

    def generate_tree_factory(context):
        def generate_tree(compsrc):
            tree = parse(compsrc)

            tree.scope.use_utility_code(
                MemoryView.get_view_utility_code(options.shared_utility_qualified_name))

            tree.scope.use_utility_code(MemoryView._get_memviewslice_declare_code())
            tree.scope.use_utility_code(MemoryView._get_typeinfo_to_format_code())
            context.include_directories.append(Code.get_utility_dir())
            return tree

        return generate_tree

    orig_cimport_from_pyx = Options.cimport_from_pyx

    def set_cimport_from_pyx(cimport_from_pyx):
        def inner(node):
            Options.cimport_from_pyx = cimport_from_pyx
            return node
        return inner

    return [
        # "cimport_from_pyx=True" to force generating __Pyx_ExportFunction
        set_cimport_from_pyx(True),
        generate_tree_factory(context),
        *Pipeline.create_pipeline(context, 'pyx', exclude_classes=()),
        Pipeline.inject_pxd_code_stage_factory(context),
        Pipeline.inject_utility_code_stage_factory(context, internalise_c_class_entries=False),
        Pipeline.inject_utility_pxd_code_stage_factory(context),
        Pipeline.abort_on_errors,
        Pipeline.generate_pyx_code_stage_factory(options, result),
        set_cimport_from_pyx(orig_cimport_from_pyx),
    ]


def generate_shared_module(options):
    Errors.init_thread()
    Errors.open_listing_file(None)

    dest_c_file = options.shared_c_file_path
    module_name = os.path.splitext(os.path.basename(dest_c_file))[0]

    context = Main.Context.from_options(options)
    scope = Symtab.ModuleScope('MemoryView', parent_module = None, context = context, is_package=False)

    with tempfile.TemporaryDirectory() as tmpdirname:
        pyx_file = os.path.join(tmpdirname, f'{module_name}.pyx')
        c_file = os.path.join(tmpdirname, f'{module_name}.c')
        with open(pyx_file, 'w'):
            pass
        source_desc = FileSourceDescriptor(pyx_file)
        comp_src = Main.CompilationSource(source_desc, EncodedString(module_name), os.getcwd())
        result = Main.create_default_resultobj(comp_src, options)

        pipeline = create_shared_library_pipeline(context, scope, options, result)
        err, enddata = Pipeline.run_pipeline(pipeline, comp_src)
        if err is None:
            shutil.copy(c_file, dest_c_file)

    return err, enddata
