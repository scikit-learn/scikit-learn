"""Calculate some properties of classes.

These happen after semantic analysis and before type checking.
"""

from typing import List, Set, Optional
from typing_extensions import Final

from mypy.nodes import (
    Node, TypeInfo, Var, Decorator, OverloadedFuncDef, SymbolTable, CallExpr, PromoteExpr,
)
from mypy.types import Instance, Type
from mypy.errors import Errors
from mypy.options import Options

# Hard coded type promotions (shared between all Python versions).
# These add extra ad-hoc edges to the subtyping relation. For example,
# int is considered a subtype of float, even though there is no
# subclass relationship.
TYPE_PROMOTIONS: Final = {
    'builtins.int': 'float',
    'builtins.float': 'complex',
}

# Hard coded type promotions for Python 3.
#
# Note that the bytearray -> bytes promotion is a little unsafe
# as some functions only accept bytes objects. Here convenience
# trumps safety.
TYPE_PROMOTIONS_PYTHON3: Final = TYPE_PROMOTIONS.copy()
TYPE_PROMOTIONS_PYTHON3.update({
    'builtins.bytearray': 'bytes',
    'builtins.memoryview': 'bytes',
})

# Hard coded type promotions for Python 2.
#
# These promotions are unsafe, but we are doing them anyway
# for convenience and also for Python 3 compatibility
# (bytearray -> str).
TYPE_PROMOTIONS_PYTHON2: Final = TYPE_PROMOTIONS.copy()
TYPE_PROMOTIONS_PYTHON2.update({
    'builtins.str': 'unicode',
    'builtins.bytearray': 'str',
    'builtins.memoryview': 'str',
})


def calculate_class_abstract_status(typ: TypeInfo, is_stub_file: bool, errors: Errors) -> None:
    """Calculate abstract status of a class.

    Set is_abstract of the type to True if the type has an unimplemented
    abstract attribute.  Also compute a list of abstract attributes.
    Report error is required ABCMeta metaclass is missing.
    """
    if typ.typeddict_type:
        return  # TypedDict can't be abstract
    concrete: Set[str] = set()
    abstract: List[str] = []
    abstract_in_this_class: List[str] = []
    if typ.is_newtype:
        # Special case: NewTypes are considered as always non-abstract, so they can be used as:
        #     Config = NewType('Config', Mapping[str, str])
        #     default = Config({'cannot': 'modify'})  # OK
        typ.abstract_attributes = []
        return
    for base in typ.mro:
        for name, symnode in base.names.items():
            node = symnode.node
            if isinstance(node, OverloadedFuncDef):
                # Unwrap an overloaded function definition. We can just
                # check arbitrarily the first overload item. If the
                # different items have a different abstract status, there
                # should be an error reported elsewhere.
                if node.items:  # can be empty for invalid overloads
                    func: Optional[Node] = node.items[0]
                else:
                    func = None
            else:
                func = node
            if isinstance(func, Decorator):
                fdef = func.func
                if fdef.is_abstract and name not in concrete:
                    typ.is_abstract = True
                    abstract.append(name)
                    if base is typ:
                        abstract_in_this_class.append(name)
            elif isinstance(node, Var):
                if node.is_abstract_var and name not in concrete:
                    typ.is_abstract = True
                    abstract.append(name)
                    if base is typ:
                        abstract_in_this_class.append(name)
            concrete.add(name)
    # In stubs, abstract classes need to be explicitly marked because it is too
    # easy to accidentally leave a concrete class abstract by forgetting to
    # implement some methods.
    typ.abstract_attributes = sorted(abstract)
    if is_stub_file:
        if typ.declared_metaclass and typ.declared_metaclass.type.fullname == 'abc.ABCMeta':
            return
        if typ.is_protocol:
            return
        if abstract and not abstract_in_this_class:
            def report(message: str, severity: str) -> None:
                errors.report(typ.line, typ.column, message, severity=severity)

            attrs = ", ".join('"{}"'.format(attr) for attr in sorted(abstract))
            report("Class {} has abstract attributes {}".format(typ.fullname, attrs), 'error')
            report("If it is meant to be abstract, add 'abc.ABCMeta' as an explicit metaclass",
                   'note')
    if typ.is_final and abstract:
        attrs = ", ".join('"{}"'.format(attr) for attr in sorted(abstract))
        errors.report(typ.line, typ.column,
                      "Final class {} has abstract attributes {}".format(typ.fullname, attrs))


def check_protocol_status(info: TypeInfo, errors: Errors) -> None:
    """Check that all classes in MRO of a protocol are protocols"""
    if info.is_protocol:
        for type in info.bases:
            if not type.type.is_protocol and type.type.fullname != 'builtins.object':
                def report(message: str, severity: str) -> None:
                    errors.report(info.line, info.column, message, severity=severity)
                report('All bases of a protocol must be protocols', 'error')


def calculate_class_vars(info: TypeInfo) -> None:
    """Try to infer additional class variables.

    Subclass attribute assignments with no type annotation are assumed
    to be classvar if overriding a declared classvar from the base
    class.

    This must happen after the main semantic analysis pass, since
    this depends on base class bodies having been fully analyzed.
    """
    for name, sym in info.names.items():
        node = sym.node
        if isinstance(node, Var) and node.info and node.is_inferred and not node.is_classvar:
            for base in info.mro[1:]:
                member = base.names.get(name)
                if (member is not None
                        and isinstance(member.node, Var)
                        and member.node.is_classvar):
                    node.is_classvar = True


def add_type_promotion(info: TypeInfo, module_names: SymbolTable, options: Options) -> None:
    """Setup extra, ad-hoc subtyping relationships between classes (promotion).

    This includes things like 'int' being compatible with 'float'.
    """
    defn = info.defn
    promote_target: Optional[Type] = None
    for decorator in defn.decorators:
        if isinstance(decorator, CallExpr):
            analyzed = decorator.analyzed
            if isinstance(analyzed, PromoteExpr):
                # _promote class decorator (undocumented feature).
                promote_target = analyzed.type
    if not promote_target:
        promotions = (TYPE_PROMOTIONS_PYTHON3 if options.python_version[0] >= 3
                      else TYPE_PROMOTIONS_PYTHON2)
        if defn.fullname in promotions:
            target_sym = module_names.get(promotions[defn.fullname])
            # With test stubs, the target may not exist.
            if target_sym:
                target_info = target_sym.node
                assert isinstance(target_info, TypeInfo)
                promote_target = Instance(target_info, [])
    defn.info._promote = promote_target
