"""Tests for stubs.

Verify that various things in stubs are consistent with how things behave at runtime.

"""

import argparse
import copy
import enum
import importlib
import inspect
import re
import sys
import types
import warnings
from functools import singledispatch
from pathlib import Path
from typing import Any, Dict, Generic, Iterator, List, Optional, Tuple, TypeVar, Union, cast

from typing_extensions import Type

import mypy.build
import mypy.modulefinder
import mypy.types
from mypy import nodes
from mypy.config_parser import parse_config_file
from mypy.options import Options
from mypy.util import FancyFormatter, bytes_to_human_readable_repr


class Missing:
    """Marker object for things that are missing (from a stub or the runtime)."""

    def __repr__(self) -> str:
        return "MISSING"


MISSING = Missing()

T = TypeVar("T")
MaybeMissing = Union[T, Missing]

_formatter = FancyFormatter(sys.stdout, sys.stderr, False)


def _style(message: str, **kwargs: Any) -> str:
    """Wrapper around mypy.util for fancy formatting."""
    kwargs.setdefault("color", "none")
    return _formatter.style(message, **kwargs)


class Error:
    def __init__(
        self,
        object_path: List[str],
        message: str,
        stub_object: MaybeMissing[nodes.Node],
        runtime_object: MaybeMissing[Any],
        *,
        stub_desc: Optional[str] = None,
        runtime_desc: Optional[str] = None
    ) -> None:
        """Represents an error found by stubtest.

        :param object_path: Location of the object with the error,
            e.g. ``["module", "Class", "method"]``
        :param message: Error message
        :param stub_object: The mypy node representing the stub
        :param runtime_object: Actual object obtained from the runtime
        :param stub_desc: Specialised description for the stub object, should you wish
        :param runtime_desc: Specialised description for the runtime object, should you wish

        """
        self.object_desc = ".".join(object_path)
        self.message = message
        self.stub_object = stub_object
        self.runtime_object = runtime_object
        self.stub_desc = stub_desc or str(getattr(stub_object, "type", stub_object))
        self.runtime_desc = runtime_desc or str(runtime_object)

    def is_missing_stub(self) -> bool:
        """Whether or not the error is for something missing from the stub."""
        return isinstance(self.stub_object, Missing)

    def is_positional_only_related(self) -> bool:
        """Whether or not the error is for something being (or not being) positional-only."""
        # TODO: This is hacky, use error codes or something more resilient
        return "leading double underscore" in self.message

    def get_description(self, concise: bool = False) -> str:
        """Returns a description of the error.

        :param concise: Whether to return a concise, one-line description

        """
        if concise:
            return _style(self.object_desc, bold=True) + " " + self.message

        stub_line = None
        stub_file: None = None
        if not isinstance(self.stub_object, Missing):
            stub_line = self.stub_object.line
        # TODO: Find a way of getting the stub file

        stub_loc_str = ""
        if stub_line:
            stub_loc_str += " at line {}".format(stub_line)
        if stub_file:
            stub_loc_str += " in file {}".format(Path(stub_file))

        runtime_line = None
        runtime_file = None
        if not isinstance(self.runtime_object, Missing):
            try:
                runtime_line = inspect.getsourcelines(self.runtime_object)[1]
            except (OSError, TypeError):
                pass
            try:
                runtime_file = inspect.getsourcefile(self.runtime_object)
            except TypeError:
                pass

        runtime_loc_str = ""
        if runtime_line:
            runtime_loc_str += " at line {}".format(runtime_line)
        if runtime_file:
            runtime_loc_str += " in file {}".format(Path(runtime_file))

        output = [
            _style("error: ", color="red", bold=True),
            _style(self.object_desc, bold=True),
            " ",
            self.message,
            "\n",
            "Stub:",
            _style(stub_loc_str, dim=True),
            "\n",
            _style(self.stub_desc + "\n", color="blue", dim=True),
            "Runtime:",
            _style(runtime_loc_str, dim=True),
            "\n",
            _style(self.runtime_desc + "\n", color="blue", dim=True),
        ]
        return "".join(output)


def test_module(module_name: str) -> Iterator[Error]:
    """Tests a given module's stub against introspecting it at runtime.

    Requires the stub to have been built already, accomplished by a call to ``build_stubs``.

    :param module_name: The module to test

    """
    stub = get_stub(module_name)
    if stub is None:
        yield Error([module_name], "failed to find stubs", MISSING, None)
        return

    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            runtime = importlib.import_module(module_name)
            # Also run the equivalent of `from module import *`
            # This could have the additional effect of loading not-yet-loaded submodules
            # mentioned in __all__
            __import__(module_name, fromlist=["*"])
    except Exception as e:
        yield Error([module_name], "failed to import: {}".format(e), stub, MISSING)
        return

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        yield from verify(stub, runtime, [module_name])


@singledispatch
def verify(
    stub: MaybeMissing[nodes.Node], runtime: MaybeMissing[Any], object_path: List[str]
) -> Iterator[Error]:
    """Entry point for comparing a stub to a runtime object.

    We use single dispatch based on the type of ``stub``.

    :param stub: The mypy node representing a part of the stub
    :param runtime: The runtime object corresponding to ``stub``

    """
    yield Error(object_path, "is an unknown mypy node", stub, runtime)


@verify.register(nodes.MypyFile)
def verify_mypyfile(
    stub: nodes.MypyFile, runtime: MaybeMissing[types.ModuleType], object_path: List[str]
) -> Iterator[Error]:
    if isinstance(runtime, Missing):
        yield Error(object_path, "is not present at runtime", stub, runtime)
        return
    if not isinstance(runtime, types.ModuleType):
        yield Error(object_path, "is not a module", stub, runtime)
        return

    # Check things in the stub
    to_check = set(
        m
        for m, o in stub.names.items()
        if not o.module_hidden and (not m.startswith("_") or hasattr(runtime, m))
    )

    def _belongs_to_runtime(r: types.ModuleType, attr: str) -> bool:
        obj = getattr(r, attr)
        obj_mod = getattr(obj, "__module__", None)
        if obj_mod is not None:
            return obj_mod == r.__name__
        return not isinstance(obj, types.ModuleType)

    runtime_public_contents = (
        runtime.__all__
        if hasattr(runtime, "__all__")
        else [
            m
            for m in dir(runtime)
            if not m.startswith("_")
            # Ensure that the object's module is `runtime`, since in the absence of __all__ we
            # don't have a good way to detect re-exports at runtime.
            and _belongs_to_runtime(runtime, m)
        ]
    )
    # Check all things declared in module's __all__, falling back to our best guess
    to_check.update(runtime_public_contents)
    to_check.difference_update({"__file__", "__doc__", "__name__", "__builtins__", "__package__"})

    for entry in sorted(to_check):
        stub_entry = stub.names[entry].node if entry in stub.names else MISSING
        if isinstance(stub_entry, nodes.MypyFile):
            # Don't recursively check exported modules, since that leads to infinite recursion
            continue
        assert stub_entry is not None
        yield from verify(
            stub_entry,
            getattr(runtime, entry, MISSING),
            object_path + [entry],
        )


@verify.register(nodes.TypeInfo)
def verify_typeinfo(
    stub: nodes.TypeInfo, runtime: MaybeMissing[Type[Any]], object_path: List[str]
) -> Iterator[Error]:
    if isinstance(runtime, Missing):
        yield Error(object_path, "is not present at runtime", stub, runtime, stub_desc=repr(stub))
        return
    if not isinstance(runtime, type):
        yield Error(object_path, "is not a type", stub, runtime, stub_desc=repr(stub))
        return

    # Check everything already defined in the stub
    to_check = set(stub.names)
    # There's a reasonable case to be made that we should always check all dunders, but it's
    # currently quite noisy. We could turn this into a denylist instead of an allowlist.
    to_check.update(
        # cast to workaround mypyc complaints
        m for m in cast(Any, vars)(runtime) if not m.startswith("_") or m in SPECIAL_DUNDERS
    )

    for entry in sorted(to_check):
        mangled_entry = entry
        if entry.startswith("__") and not entry.endswith("__"):
            mangled_entry = "_{}{}".format(stub.name, entry)
        stub_to_verify = next((t.names[entry].node for t in stub.mro if entry in t.names), MISSING)
        assert stub_to_verify is not None
        yield from verify(
            stub_to_verify,
            getattr(runtime, mangled_entry, MISSING),
            object_path + [entry],
        )


def _verify_static_class_methods(
    stub: nodes.FuncBase, runtime: Any, object_path: List[str]
) -> Iterator[str]:
    if stub.name in ("__new__", "__init_subclass__", "__class_getitem__"):
        # Special cased by Python, so don't bother checking
        return
    if inspect.isbuiltin(runtime):
        # The isinstance checks don't work reliably for builtins, e.g. datetime.datetime.now, so do
        # something a little hacky that seems to work well
        probably_class_method = isinstance(getattr(runtime, "__self__", None), type)
        if probably_class_method and not stub.is_class:
            yield "runtime is a classmethod but stub is not"
        if not probably_class_method and stub.is_class:
            yield "stub is a classmethod but runtime is not"
        return

    # Look the object up statically, to avoid binding by the descriptor protocol
    static_runtime = importlib.import_module(object_path[0])
    for entry in object_path[1:]:
        try:
            static_runtime = inspect.getattr_static(static_runtime, entry)
        except AttributeError:
            # This can happen with mangled names, ignore for now.
            # TODO: pass more information about ancestors of nodes/objects to verify, so we don't
            # have to do this hacky lookup. Would be useful in a couple other places too.
            return

    if isinstance(static_runtime, classmethod) and not stub.is_class:
        yield "runtime is a classmethod but stub is not"
    if not isinstance(static_runtime, classmethod) and stub.is_class:
        yield "stub is a classmethod but runtime is not"
    if isinstance(static_runtime, staticmethod) and not stub.is_static:
        yield "runtime is a staticmethod but stub is not"
    if not isinstance(static_runtime, staticmethod) and stub.is_static:
        yield "stub is a staticmethod but runtime is not"


def _verify_arg_name(
    stub_arg: nodes.Argument, runtime_arg: inspect.Parameter, function_name: str
) -> Iterator[str]:
    """Checks whether argument names match."""
    # Ignore exact names for most dunder methods
    if is_dunder(function_name, exclude_special=True):
        return

    def strip_prefix(s: str, prefix: str) -> str:
        return s[len(prefix):] if s.startswith(prefix) else s

    if strip_prefix(stub_arg.variable.name, "__") == runtime_arg.name:
        return

    def names_approx_match(a: str, b: str) -> bool:
        a = a.strip("_")
        b = b.strip("_")
        return a.startswith(b) or b.startswith(a) or len(a) == 1 or len(b) == 1

    # Be more permissive about names matching for positional-only arguments
    if runtime_arg.kind == inspect.Parameter.POSITIONAL_ONLY and names_approx_match(
        stub_arg.variable.name, runtime_arg.name
    ):
        return
    # This comes up with namedtuples, so ignore
    if stub_arg.variable.name == "_self":
        return
    yield (
        'stub argument "{}" differs from runtime argument "{}"'.format(
            stub_arg.variable.name, runtime_arg.name
        )
    )


def _verify_arg_default_value(
    stub_arg: nodes.Argument, runtime_arg: inspect.Parameter
) -> Iterator[str]:
    """Checks whether argument default values are compatible."""
    if runtime_arg.default != inspect.Parameter.empty:
        if stub_arg.kind.is_required():
            yield (
                'runtime argument "{}" has a default value but stub argument does not'.format(
                    runtime_arg.name
                )
            )
        else:
            runtime_type = get_mypy_type_of_runtime_value(runtime_arg.default)
            # Fallback to the type annotation type if var type is missing. The type annotation
            # is an UnboundType, but I don't know enough to know what the pros and cons here are.
            # UnboundTypes have ugly question marks following them, so default to var type.
            # Note we do this same fallback when constructing signatures in from_overloadedfuncdef
            stub_type = stub_arg.variable.type or stub_arg.type_annotation
            if isinstance(stub_type, mypy.types.TypeVarType):
                stub_type = stub_type.upper_bound
            if (
                runtime_type is not None
                and stub_type is not None
                # Avoid false positives for marker objects
                and type(runtime_arg.default) != object
                and not is_subtype_helper(runtime_type, stub_type)
            ):
                yield (
                    'runtime argument "{}" has a default value of type {}, '
                    "which is incompatible with stub argument type {}".format(
                        runtime_arg.name, runtime_type, stub_type
                    )
                )
    else:
        if stub_arg.kind.is_optional():
            yield (
                'stub argument "{}" has a default value but runtime argument does not'.format(
                    stub_arg.variable.name
                )
            )


def maybe_strip_cls(name: str, args: List[nodes.Argument]) -> List[nodes.Argument]:
    if name in ("__init_subclass__", "__class_getitem__"):
        # These are implicitly classmethods. If the stub chooses not to have @classmethod, we
        # should remove the cls argument
        if args[0].variable.name == "cls":
            return args[1:]
    return args


class Signature(Generic[T]):
    def __init__(self) -> None:
        self.pos: List[T] = []
        self.kwonly: Dict[str, T] = {}
        self.varpos: Optional[T] = None
        self.varkw: Optional[T] = None

    def __str__(self) -> str:
        def get_name(arg: Any) -> str:
            if isinstance(arg, inspect.Parameter):
                return arg.name
            if isinstance(arg, nodes.Argument):
                return arg.variable.name
            raise AssertionError

        def get_type(arg: Any) -> Optional[str]:
            if isinstance(arg, inspect.Parameter):
                return None
            if isinstance(arg, nodes.Argument):
                return str(arg.variable.type or arg.type_annotation)
            raise AssertionError

        def has_default(arg: Any) -> bool:
            if isinstance(arg, inspect.Parameter):
                return arg.default != inspect.Parameter.empty
            if isinstance(arg, nodes.Argument):
                return arg.kind.is_optional()
            raise AssertionError

        def get_desc(arg: Any) -> str:
            arg_type = get_type(arg)
            return (
                get_name(arg)
                + (": {}".format(arg_type) if arg_type else "")
                + (" = ..." if has_default(arg) else "")
            )

        kw_only = sorted(self.kwonly.values(), key=lambda a: (has_default(a), get_name(a)))
        ret = "def ("
        ret += ", ".join(
            [get_desc(arg) for arg in self.pos]
            + (["*" + get_name(self.varpos)] if self.varpos else (["*"] if self.kwonly else []))
            + [get_desc(arg) for arg in kw_only]
            + (["**" + get_name(self.varkw)] if self.varkw else [])
        )
        ret += ")"
        return ret

    @staticmethod
    def from_funcitem(stub: nodes.FuncItem) -> "Signature[nodes.Argument]":
        stub_sig: Signature[nodes.Argument] = Signature()
        stub_args = maybe_strip_cls(stub.name, stub.arguments)
        for stub_arg in stub_args:
            if stub_arg.kind.is_positional():
                stub_sig.pos.append(stub_arg)
            elif stub_arg.kind.is_named():
                stub_sig.kwonly[stub_arg.variable.name] = stub_arg
            elif stub_arg.kind == nodes.ARG_STAR:
                stub_sig.varpos = stub_arg
            elif stub_arg.kind == nodes.ARG_STAR2:
                stub_sig.varkw = stub_arg
            else:
                raise AssertionError
        return stub_sig

    @staticmethod
    def from_inspect_signature(signature: inspect.Signature) -> "Signature[inspect.Parameter]":
        runtime_sig: Signature[inspect.Parameter] = Signature()
        for runtime_arg in signature.parameters.values():
            if runtime_arg.kind in (
                inspect.Parameter.POSITIONAL_ONLY,
                inspect.Parameter.POSITIONAL_OR_KEYWORD,
            ):
                runtime_sig.pos.append(runtime_arg)
            elif runtime_arg.kind == inspect.Parameter.KEYWORD_ONLY:
                runtime_sig.kwonly[runtime_arg.name] = runtime_arg
            elif runtime_arg.kind == inspect.Parameter.VAR_POSITIONAL:
                runtime_sig.varpos = runtime_arg
            elif runtime_arg.kind == inspect.Parameter.VAR_KEYWORD:
                runtime_sig.varkw = runtime_arg
            else:
                raise AssertionError
        return runtime_sig

    @staticmethod
    def from_overloadedfuncdef(stub: nodes.OverloadedFuncDef) -> "Signature[nodes.Argument]":
        """Returns a Signature from an OverloadedFuncDef.

        If life were simple, to verify_overloadedfuncdef, we'd just verify_funcitem for each of its
        items. Unfortunately, life isn't simple and overloads are pretty deceitful. So instead, we
        try and combine the overload's items into a single signature that is compatible with any
        lies it might try to tell.

        """
        # For most dunder methods, just assume all args are positional-only
        assume_positional_only = is_dunder(stub.name, exclude_special=True)

        all_args: Dict[str, List[Tuple[nodes.Argument, int]]] = {}
        for func in map(_resolve_funcitem_from_decorator, stub.items):
            assert func is not None
            args = maybe_strip_cls(stub.name, func.arguments)
            for index, arg in enumerate(args):
                # For positional-only args, we allow overloads to have different names for the same
                # argument. To accomplish this, we just make up a fake index-based name.
                name = (
                    "__{}".format(index)
                    if arg.variable.name.startswith("__") or assume_positional_only
                    else arg.variable.name
                )
                all_args.setdefault(name, []).append((arg, index))

        def get_position(arg_name: str) -> int:
            # We just need this to return the positional args in the correct order.
            return max(index for _, index in all_args[arg_name])

        def get_type(arg_name: str) -> mypy.types.ProperType:
            with mypy.state.strict_optional_set(True):
                all_types = [
                    arg.variable.type or arg.type_annotation for arg, _ in all_args[arg_name]
                ]
                return mypy.typeops.make_simplified_union([t for t in all_types if t])

        def get_kind(arg_name: str) -> nodes.ArgKind:
            kinds = {arg.kind for arg, _ in all_args[arg_name]}
            if nodes.ARG_STAR in kinds:
                return nodes.ARG_STAR
            if nodes.ARG_STAR2 in kinds:
                return nodes.ARG_STAR2
            # The logic here is based on two tenets:
            # 1) If an arg is ever optional (or unspecified), it is optional
            # 2) If an arg is ever positional, it is positional
            is_opt = (
                len(all_args[arg_name]) < len(stub.items)
                or nodes.ARG_OPT in kinds
                or nodes.ARG_NAMED_OPT in kinds
            )
            is_pos = nodes.ARG_OPT in kinds or nodes.ARG_POS in kinds
            if is_opt:
                return nodes.ARG_OPT if is_pos else nodes.ARG_NAMED_OPT
            return nodes.ARG_POS if is_pos else nodes.ARG_NAMED

        sig: Signature[nodes.Argument] = Signature()
        for arg_name in sorted(all_args, key=get_position):
            # example_arg_name gives us a real name (in case we had a fake index-based name)
            example_arg_name = all_args[arg_name][0][0].variable.name
            arg = nodes.Argument(
                nodes.Var(example_arg_name, get_type(arg_name)),
                type_annotation=None,
                initializer=None,
                kind=get_kind(arg_name),
            )
            if arg.kind.is_positional():
                sig.pos.append(arg)
            elif arg.kind.is_named():
                sig.kwonly[arg.variable.name] = arg
            elif arg.kind == nodes.ARG_STAR:
                sig.varpos = arg
            elif arg.kind == nodes.ARG_STAR2:
                sig.varkw = arg
            else:
                raise AssertionError
        return sig


def _verify_signature(
    stub: Signature[nodes.Argument], runtime: Signature[inspect.Parameter], function_name: str
) -> Iterator[str]:
    # Check positional arguments match up
    for stub_arg, runtime_arg in zip(stub.pos, runtime.pos):
        yield from _verify_arg_name(stub_arg, runtime_arg, function_name)
        yield from _verify_arg_default_value(stub_arg, runtime_arg)
        if (
            runtime_arg.kind == inspect.Parameter.POSITIONAL_ONLY
            and not stub_arg.variable.name.startswith("__")
            and not stub_arg.variable.name.strip("_") == "self"
            and not is_dunder(function_name, exclude_special=True)  # noisy for dunder methods
        ):
            yield (
                'stub argument "{}" should be positional-only '
                '(rename with a leading double underscore, i.e. "__{}")'.format(
                    stub_arg.variable.name, runtime_arg.name
                )
            )
        if (
            runtime_arg.kind != inspect.Parameter.POSITIONAL_ONLY
            and stub_arg.variable.name.startswith("__")
        ):
            yield (
                'stub argument "{}" should be positional or keyword '
                "(remove leading double underscore)".format(stub_arg.variable.name)
            )

    # Check unmatched positional args
    if len(stub.pos) > len(runtime.pos):
        # There are cases where the stub exhaustively lists out the extra parameters the function
        # would take through *args. Hence, a) we can't check that the runtime actually takes those
        # parameters and b) below, we don't enforce that the stub takes *args, since runtime logic
        # may prevent those arguments from actually being accepted.
        if runtime.varpos is None:
            for stub_arg in stub.pos[len(runtime.pos):]:
                # If the variable is in runtime.kwonly, it's just mislabelled as not a
                # keyword-only argument
                if stub_arg.variable.name not in runtime.kwonly:
                    yield 'runtime does not have argument "{}"'.format(stub_arg.variable.name)
                else:
                    yield 'stub argument "{}" is not keyword-only'.format(stub_arg.variable.name)
            if stub.varpos is not None:
                yield 'runtime does not have *args argument "{}"'.format(stub.varpos.variable.name)
    elif len(stub.pos) < len(runtime.pos):
        for runtime_arg in runtime.pos[len(stub.pos):]:
            if runtime_arg.name not in stub.kwonly:
                yield 'stub does not have argument "{}"'.format(runtime_arg.name)
            else:
                yield 'runtime argument "{}" is not keyword-only'.format(runtime_arg.name)

    # Checks involving *args
    if len(stub.pos) <= len(runtime.pos) or runtime.varpos is None:
        if stub.varpos is None and runtime.varpos is not None:
            yield 'stub does not have *args argument "{}"'.format(runtime.varpos.name)
        if stub.varpos is not None and runtime.varpos is None:
            yield 'runtime does not have *args argument "{}"'.format(stub.varpos.variable.name)

    # Check keyword-only args
    for arg in sorted(set(stub.kwonly) & set(runtime.kwonly)):
        stub_arg, runtime_arg = stub.kwonly[arg], runtime.kwonly[arg]
        yield from _verify_arg_name(stub_arg, runtime_arg, function_name)
        yield from _verify_arg_default_value(stub_arg, runtime_arg)

    # Check unmatched keyword-only args
    if runtime.varkw is None or not set(runtime.kwonly).issubset(set(stub.kwonly)):
        # There are cases where the stub exhaustively lists out the extra parameters the function
        # would take through *kwargs. Hence, a) we only check if the runtime actually takes those
        # parameters when the above condition holds and b) below, we don't enforce that the stub
        # takes *kwargs, since runtime logic may prevent additional arguments from actually being
        # accepted.
        for arg in sorted(set(stub.kwonly) - set(runtime.kwonly)):
            yield 'runtime does not have argument "{}"'.format(arg)
    for arg in sorted(set(runtime.kwonly) - set(stub.kwonly)):
        if arg in set(stub_arg.variable.name for stub_arg in stub.pos):
            # Don't report this if we've reported it before
            if len(stub.pos) > len(runtime.pos) and runtime.varpos is not None:
                yield 'stub argument "{}" is not keyword-only'.format(arg)
        else:
            yield 'stub does not have argument "{}"'.format(arg)

    # Checks involving **kwargs
    if stub.varkw is None and runtime.varkw is not None:
        # As mentioned above, don't enforce that the stub takes **kwargs.
        # Also check against positional parameters, to avoid a nitpicky message when an argument
        # isn't marked as keyword-only
        stub_pos_names = set(stub_arg.variable.name for stub_arg in stub.pos)
        # Ideally we'd do a strict subset check, but in practice the errors from that aren't useful
        if not set(runtime.kwonly).issubset(set(stub.kwonly) | stub_pos_names):
            yield 'stub does not have **kwargs argument "{}"'.format(runtime.varkw.name)
    if stub.varkw is not None and runtime.varkw is None:
        yield 'runtime does not have **kwargs argument "{}"'.format(stub.varkw.variable.name)


@verify.register(nodes.FuncItem)
def verify_funcitem(
    stub: nodes.FuncItem, runtime: MaybeMissing[Any], object_path: List[str]
) -> Iterator[Error]:
    if isinstance(runtime, Missing):
        yield Error(object_path, "is not present at runtime", stub, runtime)
        return

    if not is_probably_a_function(runtime):
        yield Error(object_path, "is not a function", stub, runtime)
        if not callable(runtime):
            return

    for message in _verify_static_class_methods(stub, runtime, object_path):
        yield Error(object_path, "is inconsistent, " + message, stub, runtime)

    signature = safe_inspect_signature(runtime)
    if not signature:
        return

    stub_sig = Signature.from_funcitem(stub)
    runtime_sig = Signature.from_inspect_signature(signature)

    for message in _verify_signature(stub_sig, runtime_sig, function_name=stub.name):
        yield Error(
            object_path,
            "is inconsistent, " + message,
            stub,
            runtime,
            runtime_desc="def " + str(signature),
        )


@verify.register(Missing)
def verify_none(
    stub: Missing, runtime: MaybeMissing[Any], object_path: List[str]
) -> Iterator[Error]:
    yield Error(object_path, "is not present in stub", stub, runtime)


@verify.register(nodes.Var)
def verify_var(
    stub: nodes.Var, runtime: MaybeMissing[Any], object_path: List[str]
) -> Iterator[Error]:
    if isinstance(runtime, Missing):
        # Don't always yield an error here, because we often can't find instance variables
        if len(object_path) <= 2:
            yield Error(object_path, "is not present at runtime", stub, runtime)
        return

    runtime_type = get_mypy_type_of_runtime_value(runtime)
    if (
        runtime_type is not None
        and stub.type is not None
        and not is_subtype_helper(runtime_type, stub.type)
    ):
        should_error = True
        # Avoid errors when defining enums, since runtime_type is the enum itself, but we'd
        # annotate it with the type of runtime.value
        if isinstance(runtime, enum.Enum):
            runtime_type = get_mypy_type_of_runtime_value(runtime.value)
            if runtime_type is not None and is_subtype_helper(runtime_type, stub.type):
                should_error = False

        if should_error:
            yield Error(
                object_path,
                "variable differs from runtime type {}".format(runtime_type),
                stub,
                runtime,
            )


@verify.register(nodes.OverloadedFuncDef)
def verify_overloadedfuncdef(
    stub: nodes.OverloadedFuncDef, runtime: MaybeMissing[Any], object_path: List[str]
) -> Iterator[Error]:
    if isinstance(runtime, Missing):
        yield Error(object_path, "is not present at runtime", stub, runtime)
        return

    if stub.is_property:
        # We get here in cases of overloads from property.setter
        return

    if not is_probably_a_function(runtime):
        yield Error(object_path, "is not a function", stub, runtime)
        if not callable(runtime):
            return

    for message in _verify_static_class_methods(stub, runtime, object_path):
        yield Error(object_path, "is inconsistent, " + message, stub, runtime)

    signature = safe_inspect_signature(runtime)
    if not signature:
        return

    stub_sig = Signature.from_overloadedfuncdef(stub)
    runtime_sig = Signature.from_inspect_signature(signature)

    for message in _verify_signature(stub_sig, runtime_sig, function_name=stub.name):
        # TODO: This is a little hacky, but the addition here is super useful
        if "has a default value of type" in message:
            message += (
                ". This is often caused by overloads failing to account for explicitly passing "
                "in the default value."
            )
        yield Error(
            object_path,
            "is inconsistent, " + message,
            stub,
            runtime,
            stub_desc=str(stub.type) + "\nInferred signature: {}".format(stub_sig),
            runtime_desc="def " + str(signature),
        )


@verify.register(nodes.TypeVarExpr)
def verify_typevarexpr(
    stub: nodes.TypeVarExpr, runtime: MaybeMissing[Any], object_path: List[str]
) -> Iterator[Error]:
    if False:
        yield None


def _verify_property(stub: nodes.Decorator, runtime: Any) -> Iterator[str]:
    assert stub.func.is_property
    if isinstance(runtime, property):
        return
    if inspect.isdatadescriptor(runtime):
        # It's enough like a property...
        return
    # Sometimes attributes pretend to be properties, for instance, to express that they
    # are read only. So allowlist if runtime_type matches the return type of stub.
    runtime_type = get_mypy_type_of_runtime_value(runtime)
    func_type = (
        stub.func.type.ret_type if isinstance(stub.func.type, mypy.types.CallableType) else None
    )
    if (
        runtime_type is not None
        and func_type is not None
        and is_subtype_helper(runtime_type, func_type)
    ):
        return
    yield "is inconsistent, cannot reconcile @property on stub with runtime object"


def _resolve_funcitem_from_decorator(dec: nodes.OverloadPart) -> Optional[nodes.FuncItem]:
    """Returns a FuncItem that corresponds to the output of the decorator.

    Returns None if we can't figure out what that would be. For convenience, this function also
    accepts FuncItems.

    """
    if isinstance(dec, nodes.FuncItem):
        return dec
    if dec.func.is_property:
        return None

    def apply_decorator_to_funcitem(
        decorator: nodes.Expression, func: nodes.FuncItem
    ) -> Optional[nodes.FuncItem]:
        if not isinstance(decorator, nodes.RefExpr):
            return None
        if decorator.fullname is None:
            # Happens with namedtuple
            return None
        if decorator.fullname in (
            "builtins.staticmethod",
            "typing.overload",
            "abc.abstractmethod",
        ):
            return func
        if decorator.fullname == "builtins.classmethod":
            assert func.arguments[0].variable.name in ("cls", "metacls")
            ret = copy.copy(func)
            # Remove the cls argument, since it's not present in inspect.signature of classmethods
            ret.arguments = ret.arguments[1:]
            return ret
        # Just give up on any other decorators. After excluding properties, we don't run into
        # anything else when running on typeshed's stdlib.
        return None

    func: nodes.FuncItem = dec.func
    for decorator in dec.original_decorators:
        resulting_func = apply_decorator_to_funcitem(decorator, func)
        if resulting_func is None:
            return None
        func = resulting_func
    return func


@verify.register(nodes.Decorator)
def verify_decorator(
    stub: nodes.Decorator, runtime: MaybeMissing[Any], object_path: List[str]
) -> Iterator[Error]:
    if isinstance(runtime, Missing):
        yield Error(object_path, "is not present at runtime", stub, runtime)
        return
    if stub.func.is_property:
        for message in _verify_property(stub, runtime):
            yield Error(object_path, message, stub, runtime)
        return

    func = _resolve_funcitem_from_decorator(stub)
    if func is not None:
        yield from verify(func, runtime, object_path)


@verify.register(nodes.TypeAlias)
def verify_typealias(
    stub: nodes.TypeAlias, runtime: MaybeMissing[Any], object_path: List[str]
) -> Iterator[Error]:
    if isinstance(runtime, Missing):
        # ignore type aliases that don't have a runtime counterpart
        return
    stub_target = mypy.types.get_proper_type(stub.target)
    if isinstance(stub_target, mypy.types.Instance):
        yield from verify(stub_target.type, runtime, object_path)
        return
    if isinstance(stub_target, mypy.types.UnionType):
        if not getattr(runtime, "__origin__", None) is Union:
            yield Error(object_path, "is not a Union", stub, runtime, stub_desc=str(stub_target))
        # could check Union contents here...
        return
    if isinstance(stub_target, mypy.types.TupleType):
        if tuple not in getattr(runtime, "__mro__", ()):
            yield Error(
                object_path, "is not a subclass of tuple", stub, runtime,
                stub_desc=str(stub_target)
            )
        # could check Tuple contents here...
        return
    if isinstance(stub_target, mypy.types.AnyType):
        return
    yield Error(
        object_path, "is not a recognised type alias", stub, runtime, stub_desc=str(stub_target)
    )


SPECIAL_DUNDERS = ("__init__", "__new__", "__call__", "__init_subclass__", "__class_getitem__")


def is_dunder(name: str, exclude_special: bool = False) -> bool:
    """Returns whether name is a dunder name.

    :param exclude_special: Whether to return False for a couple special dunder methods.

    """
    if exclude_special and name in SPECIAL_DUNDERS:
        return False
    return name.startswith("__") and name.endswith("__")


def is_probably_a_function(runtime: Any) -> bool:
    return (
        isinstance(runtime, (types.FunctionType, types.BuiltinFunctionType))
        or isinstance(runtime, (types.MethodType, types.BuiltinMethodType))
        or (inspect.ismethoddescriptor(runtime) and callable(runtime))
    )


def safe_inspect_signature(runtime: Any) -> Optional[inspect.Signature]:
    try:
        return inspect.signature(runtime)
    except Exception:
        # inspect.signature throws ValueError all the time
        # catch RuntimeError because of https://bugs.python.org/issue39504
        # catch TypeError because of https://github.com/python/typeshed/pull/5762
        # catch AttributeError because of inspect.signature(_curses.window.border)
        return None


def is_subtype_helper(left: mypy.types.Type, right: mypy.types.Type) -> bool:
    """Checks whether ``left`` is a subtype of ``right``."""
    left = mypy.types.get_proper_type(left)
    right = mypy.types.get_proper_type(right)
    if (
        isinstance(left, mypy.types.LiteralType)
        and isinstance(left.value, int)
        and left.value in (0, 1)
        and isinstance(right, mypy.types.Instance)
        and right.type.fullname == "builtins.bool"
    ):
        # Pretend Literal[0, 1] is a subtype of bool to avoid unhelpful errors.
        return True

    with mypy.state.strict_optional_set(True):
        return mypy.subtypes.is_subtype(left, right)


def get_mypy_type_of_runtime_value(runtime: Any) -> Optional[mypy.types.Type]:
    """Returns a mypy type object representing the type of ``runtime``.

    Returns None if we can't find something that works.

    """
    if runtime is None:
        return mypy.types.NoneType()
    if isinstance(runtime, property):
        # Give up on properties to avoid issues with things that are typed as attributes.
        return None

    def anytype() -> mypy.types.AnyType:
        return mypy.types.AnyType(mypy.types.TypeOfAny.unannotated)

    if isinstance(
        runtime,
        (types.FunctionType, types.BuiltinFunctionType,
        types.MethodType, types.BuiltinMethodType)
    ):
        builtins = get_stub("builtins")
        assert builtins is not None
        type_info = builtins.names["function"].node
        assert isinstance(type_info, nodes.TypeInfo)
        fallback = mypy.types.Instance(type_info, [anytype()])
        signature = safe_inspect_signature(runtime)
        if signature:
            arg_types = []
            arg_kinds = []
            arg_names = []
            for arg in signature.parameters.values():
                arg_types.append(anytype())
                arg_names.append(
                    None if arg.kind == inspect.Parameter.POSITIONAL_ONLY else arg.name
                )
                has_default = arg.default == inspect.Parameter.empty
                if arg.kind == inspect.Parameter.POSITIONAL_ONLY:
                    arg_kinds.append(nodes.ARG_POS if has_default else nodes.ARG_OPT)
                elif arg.kind == inspect.Parameter.POSITIONAL_OR_KEYWORD:
                    arg_kinds.append(nodes.ARG_POS if has_default else nodes.ARG_OPT)
                elif arg.kind == inspect.Parameter.KEYWORD_ONLY:
                    arg_kinds.append(nodes.ARG_NAMED if has_default else nodes.ARG_NAMED_OPT)
                elif arg.kind == inspect.Parameter.VAR_POSITIONAL:
                    arg_kinds.append(nodes.ARG_STAR)
                elif arg.kind == inspect.Parameter.VAR_KEYWORD:
                    arg_kinds.append(nodes.ARG_STAR2)
                else:
                    raise AssertionError
        else:
            arg_types = [anytype(), anytype()]
            arg_kinds = [nodes.ARG_STAR, nodes.ARG_STAR2]
            arg_names = [None, None]

        return mypy.types.CallableType(
            arg_types,
            arg_kinds,
            arg_names,
            ret_type=anytype(),
            fallback=fallback,
            is_ellipsis_args=True,
        )

    # Try and look up a stub for the runtime object
    stub = get_stub(type(runtime).__module__)
    if stub is None:
        return None
    type_name = type(runtime).__name__
    if type_name not in stub.names:
        return None
    type_info = stub.names[type_name].node
    if isinstance(type_info, nodes.Var):
        return type_info.type
    if not isinstance(type_info, nodes.TypeInfo):
        return None

    if isinstance(runtime, tuple):
        # Special case tuples so we construct a valid mypy.types.TupleType
        optional_items = [get_mypy_type_of_runtime_value(v) for v in runtime]
        items = [(i if i is not None else anytype()) for i in optional_items]
        fallback = mypy.types.Instance(type_info, [anytype()])
        return mypy.types.TupleType(items, fallback)

    fallback = mypy.types.Instance(type_info, [anytype() for _ in type_info.type_vars])

    value: Union[bool, int, str]
    if isinstance(runtime, bytes):
        value = bytes_to_human_readable_repr(runtime)
    elif isinstance(runtime, enum.Enum):
        value = runtime.name
    elif isinstance(runtime, (bool, int, str)):
        value = runtime
    else:
        return fallback

    return mypy.types.LiteralType(value=value, fallback=fallback)


_all_stubs: Dict[str, nodes.MypyFile] = {}


def build_stubs(modules: List[str], options: Options, find_submodules: bool = False) -> List[str]:
    """Uses mypy to construct stub objects for the given modules.

    This sets global state that ``get_stub`` can access.

    Returns all modules we might want to check. If ``find_submodules`` is False, this is equal
    to ``modules``.

    :param modules: List of modules to build stubs for.
    :param options: Mypy options for finding and building stubs.
    :param find_submodules: Whether to attempt to find submodules of the given modules as well.

    """
    data_dir = mypy.build.default_data_dir()
    search_path = mypy.modulefinder.compute_search_paths([], options, data_dir)
    find_module_cache = mypy.modulefinder.FindModuleCache(
        search_path, fscache=None, options=options
    )

    all_modules = []
    sources = []
    for module in modules:
        all_modules.append(module)
        if not find_submodules:
            module_path = find_module_cache.find_module(module)
            if not isinstance(module_path, str):
                # test_module will yield an error later when it can't find stubs
                continue
            sources.append(mypy.modulefinder.BuildSource(module_path, module, None))
        else:
            found_sources = find_module_cache.find_modules_recursive(module)
            sources.extend(found_sources)
            all_modules.extend(s.module for s in found_sources if s.module not in all_modules)

    try:
        res = mypy.build.build(sources=sources, options=options)
    except mypy.errors.CompileError as e:
        output = [
            _style("error: ", color="red", bold=True),
            "not checking stubs due to failed mypy compile:\n",
            str(e),
        ]
        print("".join(output))
        raise RuntimeError from e
    if res.errors:
        output = [
            _style("error: ", color="red", bold=True),
            "not checking stubs due to mypy build errors:\n",
        ]
        print("".join(output) + "\n".join(res.errors))
        raise RuntimeError

    global _all_stubs
    _all_stubs = res.files

    return all_modules


def get_stub(module: str) -> Optional[nodes.MypyFile]:
    """Returns a stub object for the given module, if we've built one."""
    return _all_stubs.get(module)


def get_typeshed_stdlib_modules(custom_typeshed_dir: Optional[str]) -> List[str]:
    """Returns a list of stdlib modules in typeshed (for current Python version)."""
    stdlib_py_versions = mypy.modulefinder.load_stdlib_py_versions(custom_typeshed_dir)
    packages = set()
    # Typeshed's minimum supported Python 3 is Python 3.6
    if sys.version_info < (3, 6):
        version_info = (3, 6)
    else:
        version_info = sys.version_info[0:2]
    for module, versions in stdlib_py_versions.items():
        minver, maxver = versions
        if version_info >= minver and (maxver is None or version_info <= maxver):
            packages.add(module)

    if custom_typeshed_dir:
        typeshed_dir = Path(custom_typeshed_dir)
    else:
        typeshed_dir = Path(mypy.build.default_data_dir()) / "typeshed"
    stdlib_dir = typeshed_dir / "stdlib"

    modules = []
    for path in stdlib_dir.rglob("*.pyi"):
        if path.stem == "__init__":
            path = path.parent
        module = ".".join(path.relative_to(stdlib_dir).parts[:-1] + (path.stem,))
        if module.split(".")[0] in packages:
            modules.append(module)
    return sorted(modules)


def get_allowlist_entries(allowlist_file: str) -> Iterator[str]:
    def strip_comments(s: str) -> str:
        try:
            return s[: s.index("#")].strip()
        except ValueError:
            return s.strip()

    with open(allowlist_file) as f:
        for line in f.readlines():
            entry = strip_comments(line)
            if entry:
                yield entry


def test_stubs(args: argparse.Namespace, use_builtins_fixtures: bool = False) -> int:
    """This is stubtest! It's time to test the stubs!"""
    # Load the allowlist. This is a series of strings corresponding to Error.object_desc
    # Values in the dict will store whether we used the allowlist entry or not.
    allowlist = {
        entry: False
        for allowlist_file in args.allowlist
        for entry in get_allowlist_entries(allowlist_file)
    }
    allowlist_regexes = {entry: re.compile(entry) for entry in allowlist}

    # If we need to generate an allowlist, we store Error.object_desc for each error here.
    generated_allowlist = set()

    modules = args.modules
    if args.check_typeshed:
        assert not args.modules, "Cannot pass both --check-typeshed and a list of modules"
        modules = get_typeshed_stdlib_modules(args.custom_typeshed_dir)
        annoying_modules = {"antigravity", "this"}
        modules = [m for m in modules if m not in annoying_modules]

    assert modules, "No modules to check"

    options = Options()
    options.incremental = False
    options.custom_typeshed_dir = args.custom_typeshed_dir
    options.config_file = args.mypy_config_file
    options.use_builtins_fixtures = use_builtins_fixtures

    if options.config_file:
        def set_strict_flags() -> None:  # not needed yet
            return
        parse_config_file(options, set_strict_flags, options.config_file, sys.stdout, sys.stderr)

    try:
        modules = build_stubs(modules, options, find_submodules=not args.check_typeshed)
    except RuntimeError:
        return 1

    exit_code = 0
    for module in modules:
        for error in test_module(module):
            # Filter errors
            if args.ignore_missing_stub and error.is_missing_stub():
                continue
            if args.ignore_positional_only and error.is_positional_only_related():
                continue
            if error.object_desc in allowlist:
                allowlist[error.object_desc] = True
                continue
            is_allowlisted = False
            for w in allowlist:
                if allowlist_regexes[w].fullmatch(error.object_desc):
                    allowlist[w] = True
                    is_allowlisted = True
                    break
            if is_allowlisted:
                continue

            # We have errors, so change exit code, and output whatever necessary
            exit_code = 1
            if args.generate_allowlist:
                generated_allowlist.add(error.object_desc)
                continue
            print(error.get_description(concise=args.concise))

    # Print unused allowlist entries
    if not args.ignore_unused_allowlist:
        for w in allowlist:
            # Don't consider an entry unused if it regex-matches the empty string
            # This lets us allowlist errors that don't manifest at all on some systems
            if not allowlist[w] and not allowlist_regexes[w].fullmatch(""):
                exit_code = 1
                print("note: unused allowlist entry {}".format(w))

    # Print the generated allowlist
    if args.generate_allowlist:
        for e in sorted(generated_allowlist):
            print(e)
        exit_code = 0

    return exit_code


def parse_options(args: List[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Compares stubs to objects introspected from the runtime."
    )
    parser.add_argument("modules", nargs="*", help="Modules to test")
    parser.add_argument(
        "--concise",
        action="store_true",
        help="Makes stubtest's output more concise, one line per error",
    )
    parser.add_argument(
        "--ignore-missing-stub",
        action="store_true",
        help="Ignore errors for stub missing things that are present at runtime",
    )
    parser.add_argument(
        "--ignore-positional-only",
        action="store_true",
        help="Ignore errors for whether an argument should or shouldn't be positional-only",
    )
    parser.add_argument(
        "--allowlist",
        "--whitelist",
        action="append",
        metavar="FILE",
        default=[],
        help=(
            "Use file as an allowlist. Can be passed multiple times to combine multiple "
            "allowlists. Allowlists can be created with --generate-allowlist. Allowlists "
            "support regular expressions."
        ),
    )
    parser.add_argument(
        "--generate-allowlist",
        "--generate-whitelist",
        action="store_true",
        help="Print an allowlist (to stdout) to be used with --allowlist",
    )
    parser.add_argument(
        "--ignore-unused-allowlist",
        "--ignore-unused-whitelist",
        action="store_true",
        help="Ignore unused allowlist entries",
    )
    parser.add_argument(
        "--mypy-config-file",
        metavar="FILE",
        help=(
            "Use specified mypy config file to determine mypy plugins "
            "and mypy path"
        ),
    )
    parser.add_argument(
        "--custom-typeshed-dir", metavar="DIR", help="Use the custom typeshed in DIR"
    )
    parser.add_argument(
        "--check-typeshed", action="store_true", help="Check all stdlib modules in typeshed"
    )

    return parser.parse_args(args)


def main() -> int:
    mypy.util.check_python_version("stubtest")
    return test_stubs(parse_options(sys.argv[1:]))


if __name__ == "__main__":
    sys.exit(main())
