import inspect

from .Visitor import CythonTransform
from .StringEncoding import EncodedString
from . import Options
from . import PyrexTypes
from ..CodeWriter import ExpressionWriter
from .Errors import warning


class AnnotationWriter(ExpressionWriter):
    """
    A Cython code writer for Python expressions in argument/variable annotations.
    """
    def __init__(self, description=None):
        """description is optional. If specified it is used in
        warning messages for the nodes that don't convert to string properly.
        If not specified then no messages are generated.
        """
        ExpressionWriter.__init__(self)
        self.description = description
        self.incomplete = False

    def visit_Node(self, node):
        self.put("<???>")
        self.incomplete = True
        if self.description:
            warning(node.pos,
                    "Failed to convert code to string representation in {}".format(
                        self.description), level=1)

    def visit_LambdaNode(self, node):
        # XXX Should we do better?
        self.put("<lambda>")
        self.incomplete = True
        if self.description:
            warning(node.pos,
                    "Failed to convert lambda to string representation in {}".format(
                        self.description), level=1)

    def visit_AnnotationNode(self, node):
        self.put(node.string.value)


class EmbedSignature(CythonTransform):

    def __init__(self, context):
        super().__init__(context)
        self.class_name = None
        self.class_node = None

    def _fmt_expr(self, node):
        writer = ExpressionWriter()
        result = writer.write(node)
        # print(type(node).__name__, '-->', result)
        return result

    def _fmt_annotation(self, node):
        writer = AnnotationWriter()
        result = writer.write(node)
        # print(type(node).__name__, '-->', result)
        return result

    def _setup_format(self):
        signature_format = self.current_directives['embedsignature.format']
        self.is_format_c = signature_format == 'c'
        self.is_format_python = signature_format == 'python'
        self.is_format_clinic = signature_format == 'clinic'

    def _fmt_arg(self, arg):
        arg_doc = arg.name
        annotation = None
        defaultval = None
        if arg.is_self_arg:
            if self.is_format_clinic:
                arg_doc = '$self'
        elif arg.is_type_arg:
            if self.is_format_clinic:
                arg_doc = '$type'
        elif self.is_format_c:
            if arg.type is not PyrexTypes.py_object_type:
                arg_doc = arg.type.declaration_code(arg.name, for_display=1)
        elif self.is_format_python:
            if not arg.annotation:
                annotation = self._fmt_type(arg.type)
        if arg.annotation:
            if not self.is_format_clinic:
                annotation = self._fmt_annotation(arg.annotation)
        if arg.default:
            defaultval = self._fmt_expr(arg.default)
        if annotation:
            arg_doc = arg_doc + (': %s' % annotation)
            if defaultval:
                arg_doc = arg_doc + (' = %s' % defaultval)
        elif defaultval:
            arg_doc = arg_doc + ('=%s' % defaultval)
        return arg_doc

    def _fmt_star_arg(self, arg):
        arg_doc = arg.name
        if arg.annotation:
            if not self.is_format_clinic:
                annotation = self._fmt_annotation(arg.annotation)
                arg_doc = arg_doc + (': %s' % annotation)
        return arg_doc

    def _fmt_arglist(self, args,
                     npoargs=0, npargs=0, pargs=None,
                     nkargs=0, kargs=None,
                     hide_self=False):
        arglist = []
        for arg in args:
            if not hide_self or not arg.entry.is_self_arg:
                arg_doc = self._fmt_arg(arg)
                arglist.append(arg_doc)
        if pargs:
            arg_doc = self._fmt_star_arg(pargs)
            arglist.insert(npargs + npoargs, '*%s' % arg_doc)
        elif nkargs:
            arglist.insert(npargs + npoargs, '*')
        if npoargs:
            arglist.insert(npoargs, '/')
        if kargs:
            arg_doc = self._fmt_star_arg(kargs)
            arglist.append('**%s' % arg_doc)
        return arglist

    def _fmt_type(self, type):
        if type is PyrexTypes.py_object_type:
            return None
        elif self.is_format_c:
            code = type.declaration_code("", for_display=1)
            return code
        elif self.is_format_python:
            annotation = None
            if type.is_string:
                annotation = self.current_directives['c_string_type']
            elif type.is_numeric:
                annotation = type.py_type_name()
            if annotation is None:
                code = type.declaration_code('', for_display=1)
                annotation = code.replace(' ', '_').replace('*', 'p')
            return annotation
        return None

    def _fmt_signature(self, cls_name, func_name, args,
                       npoargs=0, npargs=0, pargs=None,
                       nkargs=0, kargs=None,
                       return_expr=None, return_type=None,
                       hide_self=False):
        arglist = self._fmt_arglist(
            args, npoargs, npargs, pargs, nkargs, kargs,
            hide_self=hide_self,
        )
        arglist_doc = ', '.join(arglist)
        func_doc = '%s(%s)' % (func_name, arglist_doc)
        if self.is_format_c and cls_name:
            func_doc = '%s.%s' % (cls_name, func_doc)
        if not self.is_format_clinic:
            ret_doc = None
            if return_expr:
                ret_doc = self._fmt_annotation(return_expr)
            elif return_type:
                ret_doc = self._fmt_type(return_type)
            if ret_doc:
                func_doc = '%s -> %s' % (func_doc, ret_doc)
        return func_doc

    def _embed_signature(self, signature, node_doc):
        if self.is_format_clinic and self.current_directives['binding']:
            return node_doc
        if node_doc:
            if self.is_format_clinic:
                docfmt = "%s\n--\n\n%s"
            else:
                docfmt = "%s\n\n%s"
            node_doc = inspect.cleandoc(node_doc)
            return docfmt % (signature, node_doc)
        else:
            if self.is_format_clinic:
                docfmt = "%s\n--\n\n"
            else:
                docfmt = "%s"
            return docfmt % signature

    def __call__(self, node):
        if not Options.docstrings:
            return node
        else:
            return super().__call__(node)

    def visit_ClassDefNode(self, node):
        oldname = self.class_name
        oldclass = self.class_node
        self.class_node = node
        try:
            # PyClassDefNode
            self.class_name = node.name
        except AttributeError:
            # CClassDefNode
            self.class_name = node.class_name
        self.visitchildren(node)
        self.class_name = oldname
        self.class_node = oldclass
        return node

    def visit_LambdaNode(self, node):
        # lambda expressions so not have signature or inner functions
        return node

    def visit_DefNode(self, node):
        if not self.current_directives['embedsignature']:
            return node
        self._setup_format()

        is_constructor = False
        hide_self = False
        if node.entry.is_special:
            is_constructor = self.class_node and node.name == '__init__'
            if is_constructor:
                class_name = None
                func_name = node.name
                if self.is_format_c:
                    func_name = self.class_name
                    hide_self = True
            else:
                class_name, func_name = self.class_name, node.name
        else:
            class_name, func_name = self.class_name, node.name

        npoargs = getattr(node, 'num_posonly_args', 0)
        nkargs = getattr(node, 'num_kwonly_args', 0)
        npargs = len(node.args) - nkargs - npoargs
        signature = self._fmt_signature(
            class_name, func_name, node.args,
            npoargs, npargs, node.star_arg,
            nkargs, node.starstar_arg,
            return_expr=node.return_type_annotation,
            return_type=None, hide_self=hide_self)
        if signature:
            if is_constructor and self.is_format_c:
                doc_holder = self.class_node.entry.type.scope
            else:
                doc_holder = node.entry
            if doc_holder.doc is not None:
                old_doc = doc_holder.doc
            elif not is_constructor and getattr(node, 'py_func', None) is not None:
                old_doc = node.py_func.entry.doc
            else:
                old_doc = None
            new_doc = self._embed_signature(signature, old_doc)
            if not node.entry.is_special or is_constructor or node.entry.wrapperbase_cname is not None:
                # TODO: the wrapperbase must be generated for __doc__ to exist;
                # however this phase is run later in the pipeline than
                # Compiler/Nodes.py:declare_pyfunction, so wrapperbase_cname
                # may already be set to None
                doc_holder.doc = EncodedString(new_doc)
            if not is_constructor and getattr(node, 'py_func', None) is not None:
                node.py_func.entry.doc = EncodedString(new_doc)
        return node

    def visit_CFuncDefNode(self, node):
        if not node.overridable:  # not cpdef FOO(...):
            return node
        if not self.current_directives['embedsignature']:
            return node
        self._setup_format()

        signature = self._fmt_signature(
            self.class_name, node.declarator.base.name,
            node.declarator.args,
            return_type=node.return_type)
        if signature:
            if node.entry.doc is not None:
                old_doc = node.entry.doc
            elif getattr(node, 'py_func', None) is not None:
                old_doc = node.py_func.entry.doc
            else:
                old_doc = None
            new_doc = self._embed_signature(signature, old_doc)
            node.entry.doc = EncodedString(new_doc)
            py_func = getattr(node, 'py_func', None)
            if py_func is not None:
                py_func.entry.doc = EncodedString(new_doc)
        return node

    def visit_PropertyNode(self, node):
        if not self.current_directives['embedsignature']:
            return node
        self._setup_format()

        entry = node.entry
        body = node.body
        prop_name = entry.name
        type_name = None
        if entry.visibility == 'public':
            if self.is_format_c:
                # property synthesised from a cdef public attribute
                type_name = entry.type.declaration_code("", for_display=1)
                if not entry.type.is_pyobject:
                    type_name = "'%s'" % type_name
                elif entry.type.is_extension_type:
                    type_name = entry.type.module_name + '.' + type_name
            elif self.is_format_python:
                type_name = self._fmt_type(entry.type)
        if type_name is None:
            for stat in body.stats:
                if stat.name != '__get__':
                    continue
                if self.is_format_c:
                    prop_name = '%s.%s' % (self.class_name, prop_name)
                ret_annotation = stat.return_type_annotation
                if ret_annotation:
                    type_name = self._fmt_annotation(ret_annotation)
        if type_name is not None :
            signature = '%s: %s' % (prop_name, type_name)
            new_doc = self._embed_signature(signature, entry.doc)
            if not self.is_format_clinic:
                entry.doc = EncodedString(new_doc)
        return node
