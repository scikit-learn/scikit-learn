import rope.base.builtins
import rope.base.codeanalyze
import rope.base.evaluate
import rope.base.libutils
import rope.base.oi.soi
import rope.base.pyscopes
from rope.base import (pynamesdef as pynames, exceptions, ast,
                       astutils, pyobjects, fscommands, arguments, utils)
from rope.base.utils import pycompat

try:
    unicode
except NameError:
    unicode = str


class PyFunction(pyobjects.PyFunction):

    def __init__(self, pycore, ast_node, parent):
        rope.base.pyobjects.AbstractFunction.__init__(self)
        rope.base.pyobjects.PyDefinedObject.__init__(
            self, pycore, ast_node, parent)
        self.arguments = self.ast_node.args
        self.parameter_pyobjects = pynames._Inferred(
            self._infer_parameters, self.get_module()._get_concluded_data())
        self.returned = pynames._Inferred(self._infer_returned)
        self.parameter_pynames = None

    def _create_structural_attributes(self):
        return {}

    def _create_concluded_attributes(self):
        return {}

    def _create_scope(self):
        return rope.base.pyscopes.FunctionScope(self.pycore, self,
                                                _FunctionVisitor)

    def _infer_parameters(self):
        pyobjects = rope.base.oi.soi.infer_parameter_objects(self)
        self._handle_special_args(pyobjects)
        return pyobjects

    def _infer_returned(self, args=None):
        return rope.base.oi.soi.infer_returned_object(self, args)

    def _handle_special_args(self, pyobjects):
        if len(pyobjects) == len(self.arguments.args):
            if self.arguments.vararg:
                pyobjects.append(rope.base.builtins.get_list())
            if self.arguments.kwarg:
                pyobjects.append(rope.base.builtins.get_dict())

    def _set_parameter_pyobjects(self, pyobjects):
        if pyobjects is not None:
            self._handle_special_args(pyobjects)
        self.parameter_pyobjects.set(pyobjects)

    def get_parameters(self):
        if self.parameter_pynames is None:
            result = {}
            for index, name in enumerate(self.get_param_names()):
                # TODO: handle tuple parameters
                result[name] = pynames.ParameterName(self, index)
            self.parameter_pynames = result
        return self.parameter_pynames

    def get_parameter(self, index):
        if index < len(self.parameter_pyobjects.get()):
            return self.parameter_pyobjects.get()[index]

    def get_returned_object(self, args):
        return self.returned.get(args)

    def get_name(self):
        return self.get_ast().name

    def get_param_names(self, special_args=True):
        # TODO: handle tuple parameters
        result = [pycompat.get_ast_arg_arg(node) for node in self.arguments.args
                  if isinstance(node, pycompat.ast_arg_type)]
        if special_args:
            if self.arguments.vararg:
                result.append(pycompat.get_ast_arg_arg(self.arguments.vararg))
            if self.arguments.kwarg:
                result.append(pycompat.get_ast_arg_arg(self.arguments.kwarg))
        return result

    def get_kind(self):
        """Get function type

        It returns one of 'function', 'method', 'staticmethod' or
        'classmethod' strs.

        """
        scope = self.parent.get_scope()
        if isinstance(self.parent, PyClass):
            for decorator in self.decorators:
                pyname = rope.base.evaluate.eval_node(scope, decorator)
                if pyname == rope.base.builtins.builtins['staticmethod']:
                    return 'staticmethod'
                if pyname == rope.base.builtins.builtins['classmethod']:
                    return 'classmethod'
            return 'method'
        return 'function'

    @property
    def decorators(self):
        try:
            return getattr(self.ast_node, 'decorator_list')
        except AttributeError:
            return getattr(self.ast_node, 'decorators', None)


class PyClass(pyobjects.PyClass):

    def __init__(self, pycore, ast_node, parent):
        self.visitor_class = _ClassVisitor
        rope.base.pyobjects.AbstractClass.__init__(self)
        rope.base.pyobjects.PyDefinedObject.__init__(
            self, pycore, ast_node, parent)
        self.parent = parent
        self._superclasses = self.get_module()._get_concluded_data()

    def get_superclasses(self):
        if self._superclasses.get() is None:
            self._superclasses.set(self._get_bases())
        return self._superclasses.get()

    def get_name(self):
        return self.get_ast().name

    def _create_concluded_attributes(self):
        result = {}
        for base in reversed(self.get_superclasses()):
            result.update(base.get_attributes())
        return result

    def _get_bases(self):
        result = []
        for base_name in self.ast_node.bases:
            base = rope.base.evaluate.eval_node(self.parent.get_scope(),
                                                base_name)
            if base is not None and \
                base.get_object().get_type() == \
                    rope.base.pyobjects.get_base_type('Type'):
                    result.append(base.get_object())
        return result

    def _create_scope(self):
        return rope.base.pyscopes.ClassScope(self.pycore, self)


class PyModule(pyobjects.PyModule):

    def __init__(self, pycore, source=None,
                 resource=None, force_errors=False):
        ignore = pycore.project.prefs.get('ignore_syntax_errors', False)
        syntax_errors = force_errors or not ignore
        self.has_errors = False
        try:
            source, node = self._init_source(pycore, source, resource)
        except exceptions.ModuleSyntaxError:
            self.has_errors = True
            if syntax_errors:
                raise
            else:
                source = '\n'
                node = ast.parse('\n')
        self.source_code = source
        self.star_imports = []
        self.visitor_class = _GlobalVisitor
        self.coding = fscommands.read_str_coding(self.source_code)
        super(PyModule, self).__init__(pycore, node, resource)

    def _init_source(self, pycore, source_code, resource):
        filename = 'string'
        if resource:
            filename = resource.path
        try:
            if source_code is None:
                source_bytes = resource.read_bytes()
                source_code = fscommands.file_data_to_unicode(source_bytes)
            else:
                if isinstance(source_code, unicode):
                    source_bytes = fscommands.unicode_to_file_data(source_code)
                else:
                    source_bytes = source_code
            ast_node = ast.parse(source_bytes, filename=filename)
        except SyntaxError as e:
            raise exceptions.ModuleSyntaxError(filename, e.lineno, e.msg)
        except UnicodeDecodeError as e:
            raise exceptions.ModuleSyntaxError(filename, 1, '%s' % (e.reason))
        return source_code, ast_node

    @utils.prevent_recursion(lambda: {})
    def _create_concluded_attributes(self):
        result = {}
        for star_import in self.star_imports:
            result.update(star_import.get_names())
        return result

    def _create_scope(self):
        return rope.base.pyscopes.GlobalScope(self.pycore, self)

    @property
    @utils.saveit
    def lines(self):
        """A `SourceLinesAdapter`"""
        return rope.base.codeanalyze.SourceLinesAdapter(self.source_code)

    @property
    @utils.saveit
    def logical_lines(self):
        """A `LogicalLinesFinder`"""
        return rope.base.codeanalyze.CachingLogicalLineFinder(self.lines)

    def get_name(self):
        return rope.base.libutils.modname(self.get_resource())

class PyPackage(pyobjects.PyPackage):

    def __init__(self, pycore, resource=None, force_errors=False):
        self.resource = resource
        init_dot_py = self._get_init_dot_py()
        if init_dot_py is not None:
            ast_node = pycore.project.get_pymodule(
                init_dot_py, force_errors=force_errors).get_ast()
        else:
            ast_node = ast.parse('\n')
        super(PyPackage, self).__init__(pycore, ast_node, resource)

    def _create_structural_attributes(self):
        result = {}
        modname = rope.base.libutils.modname(self.resource)
        extension_submodules = self.pycore._builtin_submodules(modname)
        for name, module in extension_submodules.items():
            result[name] = rope.base.builtins.BuiltinName(module)
        if self.resource is None:
            return result
        for name, resource in self._get_child_resources().items():
            result[name] = pynames.ImportedModule(self, resource=resource)
        return result

    def _create_concluded_attributes(self):
        result = {}
        init_dot_py = self._get_init_dot_py()
        if init_dot_py:
            init_object = self.pycore.project.get_pymodule(init_dot_py)
            result.update(init_object.get_attributes())
        return result

    def _get_child_resources(self):
        result = {}
        for child in self.resource.get_children():
            if child.is_folder():
                result[child.name] = child
            elif child.name.endswith('.py') and \
                    child.name != '__init__.py':
                name = child.name[:-3]
                result[name] = child
        return result

    def _get_init_dot_py(self):
        if self.resource is not None and \
                self.resource.has_child('__init__.py'):
            return self.resource.get_child('__init__.py')
        else:
            return None

    def _create_scope(self):
        return self.get_module().get_scope()

    def get_module(self):
        init_dot_py = self._get_init_dot_py()
        if init_dot_py:
            return self.pycore.project.get_pymodule(init_dot_py)
        return self


class _AssignVisitor(object):

    def __init__(self, scope_visitor):
        self.scope_visitor = scope_visitor
        self.assigned_ast = None

    def _Assign(self, node):
        self.assigned_ast = node.value
        for child_node in node.targets:
            ast.walk(child_node, self)

    def _assigned(self, name, assignment=None):
        self.scope_visitor._assigned(name, assignment)

    def _Name(self, node):
        assignment = None
        if self.assigned_ast is not None:
            assignment = pynames.AssignmentValue(self.assigned_ast)
        self._assigned(node.id, assignment)

    def _Tuple(self, node):
        names = astutils.get_name_levels(node)
        for name, levels in names:
            assignment = None
            if self.assigned_ast is not None:
                assignment = pynames.AssignmentValue(self.assigned_ast, levels)
            self._assigned(name, assignment)

    def _Attribute(self, node):
        pass

    def _Subscript(self, node):
        pass

    def _Slice(self, node):
        pass


class _ScopeVisitor(object):

    def __init__(self, pycore, owner_object):
        self.pycore = pycore
        self.owner_object = owner_object
        self.names = {}
        self.defineds = []

    def get_module(self):
        if self.owner_object is not None:
            return self.owner_object.get_module()
        else:
            return None

    def _ClassDef(self, node):
        pyclass = PyClass(self.pycore, node, self.owner_object)
        self.names[node.name] = pynames.DefinedName(pyclass)
        self.defineds.append(pyclass)

    def _FunctionDef(self, node):
        pyfunction = PyFunction(self.pycore, node, self.owner_object)
        for decorator in pyfunction.decorators:
            if isinstance(decorator, ast.Name) and decorator.id == 'property':
                if isinstance(self, _ClassVisitor):
                    type_ = rope.base.builtins.Property(pyfunction)
                    arg = pynames.UnboundName(
                        rope.base.pyobjects.PyObject(self.owner_object))

                    def _eval(type_=type_, arg=arg):
                        return type_.get_property_object(
                            arguments.ObjectArguments([arg]))
                    self.names[node.name] = pynames.EvaluatedName(
                        _eval, module=self.get_module(), lineno=node.lineno)
                    break
        else:
            self.names[node.name] = pynames.DefinedName(pyfunction)
        self.defineds.append(pyfunction)

    def _Assign(self, node):
        ast.walk(node, _AssignVisitor(self))

    def _AugAssign(self, node):
        pass

    def _For(self, node):
        names = self._update_evaluated(node.target, node.iter,  # noqa
                                       '.__iter__().next()')
        for child in node.body + node.orelse:
            ast.walk(child, self)

    def _assigned(self, name, assignment):
        pyname = self.names.get(name, None)
        if pyname is None:
            pyname = pynames.AssignedName(module=self.get_module())
        if isinstance(pyname, pynames.AssignedName):
            if assignment is not None:
                pyname.assignments.append(assignment)
            self.names[name] = pyname

    def _update_evaluated(self, targets, assigned,
                          evaluation='', eval_type=False):
        result = {}
        if isinstance(targets, str):
            assignment = pynames.AssignmentValue(assigned, [],
                                                 evaluation, eval_type)
            self._assigned(targets, assignment)
        else:
            names = astutils.get_name_levels(targets)
            for name, levels in names:
                assignment = pynames.AssignmentValue(assigned, levels,
                                                     evaluation, eval_type)
                self._assigned(name, assignment)
        return result

    def _With(self, node):
        for item in pycompat.get_ast_with_items(node):
            if item.optional_vars:
                self._update_evaluated(item.optional_vars,
                                       item.context_expr, '.__enter__()')
        for child in node.body:
            ast.walk(child, self)

    def _excepthandler(self, node):
        node_name_type = str if pycompat.PY3 else ast.Name
        if node.name is not None and isinstance(node.name, node_name_type):
            type_node = node.type
            if isinstance(node.type, ast.Tuple) and type_node.elts:
                type_node = type_node.elts[0]
            self._update_evaluated(node.name, type_node, eval_type=True)

        for child in node.body:
            ast.walk(child, self)

    def _ExceptHandler(self, node):
        self._excepthandler(node)

    def _Import(self, node):
        for import_pair in node.names:
            module_name = import_pair.name
            alias = import_pair.asname
            first_package = module_name.split('.')[0]
            if alias is not None:
                imported = pynames.ImportedModule(self.get_module(),
                                                  module_name)
                if not self._is_ignored_import(imported):
                    self.names[alias] = imported
            else:
                imported = pynames.ImportedModule(self.get_module(),
                                                  first_package)
                if not self._is_ignored_import(imported):
                    self.names[first_package] = imported

    def _ImportFrom(self, node):
        level = 0
        if node.level:
            level = node.level
        imported_module = pynames.ImportedModule(self.get_module(),
                                                 node.module, level)
        if self._is_ignored_import(imported_module):
            return
        if len(node.names) == 1 and node.names[0].name == '*':
            if isinstance(self.owner_object, PyModule):
                self.owner_object.star_imports.append(
                    StarImport(imported_module))
        else:
            for imported_name in node.names:
                imported = imported_name.name
                alias = imported_name.asname
                if alias is not None:
                    imported = alias
                self.names[imported] = pynames.ImportedName(imported_module,
                                                            imported_name.name)

    def _is_ignored_import(self, imported_module):
        if not self.pycore.project.prefs.get('ignore_bad_imports', False):
            return False
        return not isinstance(imported_module.get_object(),
                              rope.base.pyobjects.AbstractModule)

    def _Global(self, node):
        module = self.get_module()
        for name in node.names:
            if module is not None:
                try:
                    pyname = module[name]
                except exceptions.AttributeNotFoundError:
                    pyname = pynames.AssignedName(node.lineno)
            self.names[name] = pyname


class _GlobalVisitor(_ScopeVisitor):

    def __init__(self, pycore, owner_object):
        super(_GlobalVisitor, self).__init__(pycore, owner_object)


class _ClassVisitor(_ScopeVisitor):

    def __init__(self, pycore, owner_object):
        super(_ClassVisitor, self).__init__(pycore, owner_object)

    def _FunctionDef(self, node):
        _ScopeVisitor._FunctionDef(self, node)
        if len(node.args.args) > 0:
            first = node.args.args[0]
            new_visitor = None
            if isinstance(first, pycompat.ast_arg_type):
                new_visitor = _ClassInitVisitor(self, pycompat.get_ast_arg_arg(first))
            if new_visitor is not None:
                for child in ast.get_child_nodes(node):
                    ast.walk(child, new_visitor)


class _FunctionVisitor(_ScopeVisitor):

    def __init__(self, pycore, owner_object):
        super(_FunctionVisitor, self).__init__(pycore, owner_object)
        self.returned_asts = []
        self.generator = False

    def _Return(self, node):
        if node.value is not None:
            self.returned_asts.append(node.value)

    def _Yield(self, node):
        if node.value is not None:
            self.returned_asts.append(node.value)
        self.generator = True


class _ClassInitVisitor(_AssignVisitor):

    def __init__(self, scope_visitor, self_name):
        super(_ClassInitVisitor, self).__init__(scope_visitor)
        self.self_name = self_name

    def _Attribute(self, node):
        if not isinstance(node.ctx, ast.Store):
            return
        if isinstance(node.value, ast.Name) and \
           node.value.id == self.self_name:
            if node.attr not in self.scope_visitor.names:
                self.scope_visitor.names[node.attr] = pynames.AssignedName(
                    lineno=node.lineno, module=self.scope_visitor.get_module())
            if self.assigned_ast is not None:
                pyname = self.scope_visitor.names[node.attr]
                if isinstance(pyname, pynames.AssignedName):
                    pyname.assignments.append(
                        pynames.AssignmentValue(self.assigned_ast))

    def _Tuple(self, node):
        if not isinstance(node.ctx, ast.Store):
            return
        for child in ast.get_child_nodes(node):
            ast.walk(child, self)

    def _Name(self, node):
        pass

    def _FunctionDef(self, node):
        pass

    def _ClassDef(self, node):
        pass

    def _For(self, node):
        pass

    def _With(self, node):
        pass


class StarImport(object):

    def __init__(self, imported_module):
        self.imported_module = imported_module

    def get_names(self):
        result = {}
        imported = self.imported_module.get_object()
        for name in imported:
            if not name.startswith('_'):
                result[name] = pynames.ImportedName(self.imported_module, name)
        return result
