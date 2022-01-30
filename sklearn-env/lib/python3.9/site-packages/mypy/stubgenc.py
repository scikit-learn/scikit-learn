#!/usr/bin/env python3
"""Stub generator for C modules.

The public interface is via the mypy.stubgen module.
"""

import importlib
import inspect
import os.path
import re
from typing import List, Dict, Tuple, Optional, Mapping, Any, Set
from types import ModuleType
from typing_extensions import Final

from mypy.moduleinspect import is_c_module
from mypy.stubdoc import (
    infer_sig_from_docstring, infer_prop_type_from_docstring, ArgSig,
    infer_arg_sig_from_anon_docstring, infer_ret_type_sig_from_anon_docstring,
    infer_ret_type_sig_from_docstring, FunctionSig
)

# Members of the typing module to consider for importing by default.
_DEFAULT_TYPING_IMPORTS: Final = (
    'Any',
    'Callable',
    'ClassVar',
    'Dict',
    'Iterable',
    'Iterator',
    'List',
    'Optional',
    'Tuple',
    'Union',
)


def generate_stub_for_c_module(module_name: str,
                               target: str,
                               sigs: Optional[Dict[str, str]] = None,
                               class_sigs: Optional[Dict[str, str]] = None) -> None:
    """Generate stub for C module.

    This combines simple runtime introspection (looking for docstrings and attributes
    with simple builtin types) and signatures inferred from .rst documentation (if given).

    If directory for target doesn't exist it will be created. Existing stub
    will be overwritten.
    """
    module = importlib.import_module(module_name)
    assert is_c_module(module), '%s is not a C module' % module_name
    subdir = os.path.dirname(target)
    if subdir and not os.path.isdir(subdir):
        os.makedirs(subdir)
    imports: List[str] = []
    functions: List[str] = []
    done = set()
    items = sorted(module.__dict__.items(), key=lambda x: x[0])
    for name, obj in items:
        if is_c_function(obj):
            generate_c_function_stub(module, name, obj, functions, imports=imports, sigs=sigs)
            done.add(name)
    types: List[str] = []
    for name, obj in items:
        if name.startswith('__') and name.endswith('__'):
            continue
        if is_c_type(obj):
            generate_c_type_stub(module, name, obj, types, imports=imports, sigs=sigs,
                                 class_sigs=class_sigs)
            done.add(name)
    variables = []
    for name, obj in items:
        if name.startswith('__') and name.endswith('__'):
            continue
        if name not in done and not inspect.ismodule(obj):
            type_str = strip_or_import(get_type_fullname(type(obj)), module, imports)
            variables.append('%s: %s' % (name, type_str))
    output = []
    for line in sorted(set(imports)):
        output.append(line)
    for line in variables:
        output.append(line)
    for line in types:
        if line.startswith('class') and output and output[-1]:
            output.append('')
        output.append(line)
    if output and functions:
        output.append('')
    for line in functions:
        output.append(line)
    output = add_typing_import(output)
    with open(target, 'w') as file:
        for line in output:
            file.write('%s\n' % line)


def add_typing_import(output: List[str]) -> List[str]:
    """Add typing imports for collections/types that occur in the generated stub."""
    names = []
    for name in _DEFAULT_TYPING_IMPORTS:
        if any(re.search(r'\b%s\b' % name, line) for line in output):
            names.append(name)
    if names:
        return ['from typing import %s' % ', '.join(names), ''] + output
    else:
        return output[:]


def is_c_function(obj: object) -> bool:
    return inspect.isbuiltin(obj) or type(obj) is type(ord)


def is_c_method(obj: object) -> bool:
    return inspect.ismethoddescriptor(obj) or type(obj) in (type(str.index),
                                                            type(str.__add__),
                                                            type(str.__new__))


def is_c_classmethod(obj: object) -> bool:
    return inspect.isbuiltin(obj) or type(obj).__name__ in ('classmethod',
                                                            'classmethod_descriptor')


def is_c_property(obj: object) -> bool:
    return inspect.isdatadescriptor(obj) and hasattr(obj, 'fget')


def is_c_property_readonly(prop: Any) -> bool:
    return prop.fset is None


def is_c_type(obj: object) -> bool:
    return inspect.isclass(obj) or type(obj) is type(int)


def is_pybind11_overloaded_function_docstring(docstr: str, name: str) -> bool:
    return docstr.startswith("{}(*args, **kwargs)\n".format(name) +
                             "Overloaded function.\n\n")


def generate_c_function_stub(module: ModuleType,
                             name: str,
                             obj: object,
                             output: List[str],
                             imports: List[str],
                             self_var: Optional[str] = None,
                             sigs: Optional[Dict[str, str]] = None,
                             class_name: Optional[str] = None,
                             class_sigs: Optional[Dict[str, str]] = None) -> None:
    """Generate stub for a single function or method.

    The result (always a single line) will be appended to 'output'.
    If necessary, any required names will be added to 'imports'.
    The 'class_name' is used to find signature of __init__ or __new__ in
    'class_sigs'.
    """
    if sigs is None:
        sigs = {}
    if class_sigs is None:
        class_sigs = {}

    ret_type = 'None' if name == '__init__' and class_name else 'Any'

    if (
        name in ("__new__", "__init__")
        and name not in sigs
        and class_name
        and class_name in class_sigs
    ):
        inferred: Optional[List[FunctionSig]] = [
            FunctionSig(
                name=name,
                args=infer_arg_sig_from_anon_docstring(class_sigs[class_name]),
                ret_type=ret_type,
            )
        ]
    else:
        docstr = getattr(obj, '__doc__', None)
        inferred = infer_sig_from_docstring(docstr, name)
        if inferred:
            assert docstr is not None
            if is_pybind11_overloaded_function_docstring(docstr, name):
                # Remove pybind11 umbrella (*args, **kwargs) for overloaded functions
                del inferred[-1]
        if not inferred:
            if class_name and name not in sigs:
                inferred = [FunctionSig(name, args=infer_method_sig(name, self_var),
                                        ret_type=ret_type)]
            else:
                inferred = [FunctionSig(name=name,
                                        args=infer_arg_sig_from_anon_docstring(
                                            sigs.get(name, '(*args, **kwargs)')),
                                        ret_type=ret_type)]
        elif class_name and self_var:
            args = inferred[0].args
            if not args or args[0].name != self_var:
                args.insert(0, ArgSig(name=self_var))

    is_overloaded = len(inferred) > 1 if inferred else False
    if is_overloaded:
        imports.append('from typing import overload')
    if inferred:
        for signature in inferred:
            sig = []
            for arg in signature.args:
                if arg.name == self_var:
                    arg_def = self_var
                else:
                    arg_def = arg.name
                    if arg_def == 'None':
                        arg_def = '_none'  # None is not a valid argument name

                    if arg.type:
                        arg_def += ": " + strip_or_import(arg.type, module, imports)

                    if arg.default:
                        arg_def += " = ..."

                sig.append(arg_def)

            if is_overloaded:
                output.append('@overload')
            output.append('def {function}({args}) -> {ret}: ...'.format(
                function=name,
                args=", ".join(sig),
                ret=strip_or_import(signature.ret_type, module, imports)
            ))


def strip_or_import(typ: str, module: ModuleType, imports: List[str]) -> str:
    """Strips unnecessary module names from typ.

    If typ represents a type that is inside module or is a type coming from builtins, remove
    module declaration from it. Return stripped name of the type.

    Arguments:
        typ: name of the type
        module: in which this type is used
        imports: list of import statements (may be modified during the call)
    """
    stripped_type = typ
    if any(c in typ for c in '[,'):
        for subtyp in re.split(r'[\[,\]]', typ):
            strip_or_import(subtyp.strip(), module, imports)
        if module:
            stripped_type = re.sub(
                r'(^|[\[, ]+)' + re.escape(module.__name__ + '.'),
                r'\1',
                typ,
            )
    elif module and typ.startswith(module.__name__ + '.'):
        stripped_type = typ[len(module.__name__) + 1:]
    elif '.' in typ:
        arg_module = typ[:typ.rindex('.')]
        if arg_module == 'builtins':
            stripped_type = typ[len('builtins') + 1:]
        else:
            imports.append('import %s' % (arg_module,))
    if stripped_type == 'NoneType':
        stripped_type = 'None'
    return stripped_type


def is_static_property(obj: object) -> bool:
    return type(obj).__name__ == 'pybind11_static_property'


def generate_c_property_stub(name: str, obj: object,
                             static_properties: List[str],
                             rw_properties: List[str],
                             ro_properties: List[str], readonly: bool,
                             module: Optional[ModuleType] = None,
                             imports: Optional[List[str]] = None) -> None:
    """Generate property stub using introspection of 'obj'.

    Try to infer type from docstring, append resulting lines to 'output'.
    """

    def infer_prop_type(docstr: Optional[str]) -> Optional[str]:
        """Infer property type from docstring or docstring signature."""
        if docstr is not None:
            inferred = infer_ret_type_sig_from_anon_docstring(docstr)
            if not inferred:
                inferred = infer_ret_type_sig_from_docstring(docstr, name)
            if not inferred:
                inferred = infer_prop_type_from_docstring(docstr)
            return inferred
        else:
            return None

    inferred = infer_prop_type(getattr(obj, '__doc__', None))
    if not inferred:
        fget = getattr(obj, 'fget', None)
        inferred = infer_prop_type(getattr(fget, '__doc__', None))
    if not inferred:
        inferred = 'Any'

    if module is not None and imports is not None:
        inferred = strip_or_import(inferred, module, imports)

    if is_static_property(obj):
        trailing_comment = "  # read-only" if readonly else ""
        static_properties.append(
            '{}: ClassVar[{}] = ...{}'.format(name, inferred, trailing_comment)
        )
    else:  # regular property
        if readonly:
            ro_properties.append('@property')
            ro_properties.append('def {}(self) -> {}: ...'.format(name, inferred))
        else:
            rw_properties.append('{}: {}'.format(name, inferred))


def generate_c_type_stub(module: ModuleType,
                         class_name: str,
                         obj: type,
                         output: List[str],
                         imports: List[str],
                         sigs: Optional[Dict[str, str]] = None,
                         class_sigs: Optional[Dict[str, str]] = None) -> None:
    """Generate stub for a single class using runtime introspection.

    The result lines will be appended to 'output'. If necessary, any
    required names will be added to 'imports'.
    """
    # typeshed gives obj.__dict__ the not quite correct type Dict[str, Any]
    # (it could be a mappingproxy!), which makes mypyc mad, so obfuscate it.
    obj_dict: Mapping[str, Any] = getattr(obj, "__dict__")  # noqa
    items = sorted(obj_dict.items(), key=lambda x: method_name_sort_key(x[0]))
    methods: List[str] = []
    types: List[str] = []
    static_properties: List[str] = []
    rw_properties: List[str] = []
    ro_properties: List[str] = []
    done: Set[str] = set()
    for attr, value in items:
        if is_c_method(value) or is_c_classmethod(value):
            done.add(attr)
            if not is_skipped_attribute(attr):
                if attr == '__new__':
                    # TODO: We should support __new__.
                    if '__init__' in obj_dict:
                        # Avoid duplicate functions if both are present.
                        # But is there any case where .__new__() has a
                        # better signature than __init__() ?
                        continue
                    attr = '__init__'
                if is_c_classmethod(value):
                    methods.append('@classmethod')
                    self_var = 'cls'
                else:
                    self_var = 'self'
                generate_c_function_stub(module, attr, value, methods, imports=imports,
                                         self_var=self_var, sigs=sigs, class_name=class_name,
                                         class_sigs=class_sigs)
        elif is_c_property(value):
            done.add(attr)
            generate_c_property_stub(attr, value, static_properties, rw_properties, ro_properties,
                                     is_c_property_readonly(value),
                                     module=module, imports=imports)
        elif is_c_type(value):
            generate_c_type_stub(module, attr, value, types, imports=imports, sigs=sigs,
                                 class_sigs=class_sigs)
            done.add(attr)

    for attr, value in items:
        if is_skipped_attribute(attr):
            continue
        if attr not in done:
            static_properties.append('%s: ClassVar[%s] = ...' % (
                attr, strip_or_import(get_type_fullname(type(value)), module, imports)))
    all_bases = type.mro(obj)
    if all_bases[-1] is object:
        # TODO: Is this always object?
        del all_bases[-1]
    # remove pybind11_object. All classes generated by pybind11 have pybind11_object in their MRO,
    # which only overrides a few functions in object type
    if all_bases and all_bases[-1].__name__ == 'pybind11_object':
        del all_bases[-1]
    # remove the class itself
    all_bases = all_bases[1:]
    # Remove base classes of other bases as redundant.
    bases: List[type] = []
    for base in all_bases:
        if not any(issubclass(b, base) for b in bases):
            bases.append(base)
    if bases:
        bases_str = '(%s)' % ', '.join(
            strip_or_import(
                get_type_fullname(base),
                module,
                imports
            ) for base in bases
        )
    else:
        bases_str = ''
    if types or static_properties or rw_properties or methods or ro_properties:
        output.append('class %s%s:' % (class_name, bases_str))
        for line in types:
            if output and output[-1] and \
                    not output[-1].startswith('class') and line.startswith('class'):
                output.append('')
            output.append('    ' + line)
        for line in static_properties:
            output.append('    %s' % line)
        for line in rw_properties:
            output.append('    %s' % line)
        for line in methods:
            output.append('    %s' % line)
        for line in ro_properties:
            output.append('    %s' % line)
    else:
        output.append('class %s%s: ...' % (class_name, bases_str))


def get_type_fullname(typ: type) -> str:
    return '%s.%s' % (typ.__module__, getattr(typ, '__qualname__', typ.__name__))


def method_name_sort_key(name: str) -> Tuple[int, str]:
    """Sort methods in classes in a typical order.

    I.e.: constructor, normal methods, special methods.
    """
    if name in ('__new__', '__init__'):
        return 0, name
    if name.startswith('__') and name.endswith('__'):
        return 2, name
    return 1, name


def is_pybind_skipped_attribute(attr: str) -> bool:
    return attr.startswith("__pybind11_module_local_")


def is_skipped_attribute(attr: str) -> bool:
    return (attr in ('__getattribute__',
                     '__str__',
                     '__repr__',
                     '__doc__',
                     '__dict__',
                     '__module__',
                     '__weakref__')  # For pickling
            or is_pybind_skipped_attribute(attr)
            )


def infer_method_sig(name: str, self_var: Optional[str] = None) -> List[ArgSig]:
    args: Optional[List[ArgSig]] = None
    if name.startswith('__') and name.endswith('__'):
        name = name[2:-2]
        if name in ('hash', 'iter', 'next', 'sizeof', 'copy', 'deepcopy', 'reduce', 'getinitargs',
                    'int', 'float', 'trunc', 'complex', 'bool', 'abs', 'bytes', 'dir', 'len',
                    'reversed', 'round', 'index', 'enter'):
            args = []
        elif name == 'getitem':
            args = [ArgSig(name='index')]
        elif name == 'setitem':
            args = [ArgSig(name='index'),
                    ArgSig(name='object')]
        elif name in ('delattr', 'getattr'):
            args = [ArgSig(name='name')]
        elif name == 'setattr':
            args = [ArgSig(name='name'),
                    ArgSig(name='value')]
        elif name == 'getstate':
            args = []
        elif name == 'setstate':
            args = [ArgSig(name='state')]
        elif name in ('eq', 'ne', 'lt', 'le', 'gt', 'ge',
                      'add', 'radd', 'sub', 'rsub', 'mul', 'rmul',
                      'mod', 'rmod', 'floordiv', 'rfloordiv', 'truediv', 'rtruediv',
                      'divmod', 'rdivmod', 'pow', 'rpow',
                      'xor', 'rxor', 'or', 'ror', 'and', 'rand', 'lshift', 'rlshift',
                      'rshift', 'rrshift',
                      'contains', 'delitem',
                      'iadd', 'iand', 'ifloordiv', 'ilshift', 'imod', 'imul', 'ior',
                      'ipow', 'irshift', 'isub', 'itruediv', 'ixor'):
            args = [ArgSig(name='other')]
        elif name in ('neg', 'pos', 'invert'):
            args = []
        elif name == 'get':
            args = [ArgSig(name='instance'),
                    ArgSig(name='owner')]
        elif name == 'set':
            args = [ArgSig(name='instance'),
                    ArgSig(name='value')]
        elif name == 'reduce_ex':
            args = [ArgSig(name='protocol')]
        elif name == 'exit':
            args = [ArgSig(name='type'),
                    ArgSig(name='value'),
                    ArgSig(name='traceback')]
    if args is None:
        args = [ArgSig(name='*args'),
                ArgSig(name='**kwargs')]
    return [ArgSig(name=self_var or 'self')] + args
