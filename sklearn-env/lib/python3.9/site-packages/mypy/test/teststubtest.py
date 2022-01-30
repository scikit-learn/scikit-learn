import contextlib
import inspect
import io
import os
import re
import sys
import tempfile
import textwrap
import unittest
from pathlib import Path
from typing import Any, Callable, Iterator, List, Optional

import mypy.stubtest
from mypy.stubtest import parse_options, test_stubs
from mypy.test.data import root_dir


@contextlib.contextmanager
def use_tmp_dir() -> Iterator[None]:
    current = os.getcwd()
    with tempfile.TemporaryDirectory() as tmp:
        try:
            os.chdir(tmp)
            yield
        finally:
            os.chdir(current)


TEST_MODULE_NAME = "test_module"

stubtest_builtins_stub = """
from typing import Generic, Mapping, Sequence, TypeVar, overload

T = TypeVar('T')
T_co = TypeVar('T_co', covariant=True)
KT = TypeVar('KT')
VT = TypeVar('VT')

class object:
    def __init__(self) -> None: pass
class type: ...

class tuple(Sequence[T_co], Generic[T_co]): ...
class dict(Mapping[KT, VT]): ...

class function: pass
class ellipsis: pass

class int: ...
class float: ...
class bool(int): ...
class str: ...
class bytes: ...

class list(Sequence[T]): ...

def property(f: T) -> T: ...
def classmethod(f: T) -> T: ...
def staticmethod(f: T) -> T: ...
"""


def run_stubtest(
    stub: str, runtime: str, options: List[str], config_file: Optional[str] = None,
) -> str:
    with use_tmp_dir():
        with open("builtins.pyi", "w") as f:
            f.write(stubtest_builtins_stub)
        with open("{}.pyi".format(TEST_MODULE_NAME), "w") as f:
            f.write(stub)
        with open("{}.py".format(TEST_MODULE_NAME), "w") as f:
            f.write(runtime)
        if config_file:
            with open("{}_config.ini".format(TEST_MODULE_NAME), "w") as f:
                f.write(config_file)
            options = options + ["--mypy-config-file", "{}_config.ini".format(TEST_MODULE_NAME)]
        if sys.path[0] != ".":
            sys.path.insert(0, ".")
        if TEST_MODULE_NAME in sys.modules:
            del sys.modules[TEST_MODULE_NAME]

        output = io.StringIO()
        with contextlib.redirect_stdout(output):
            test_stubs(
                parse_options([TEST_MODULE_NAME] + options),
                use_builtins_fixtures=True
            )

        module_path = Path(os.getcwd()) / TEST_MODULE_NAME
        # remove cwd as it's not available from outside
        return output.getvalue().replace(str(module_path), TEST_MODULE_NAME)


class Case:
    def __init__(self, stub: str, runtime: str, error: Optional[str]):
        self.stub = stub
        self.runtime = runtime
        self.error = error


def collect_cases(fn: Callable[..., Iterator[Case]]) -> Callable[..., None]:
    """run_stubtest used to be slow, so we used this decorator to combine cases.

    If you're reading this and bored, feel free to refactor this and make it more like
    other mypy tests.

    """

    def test(*args: Any, **kwargs: Any) -> None:
        cases = list(fn(*args, **kwargs))
        expected_errors = set()
        for c in cases:
            if c.error is None:
                continue
            expected_error = "{}.{}".format(TEST_MODULE_NAME, c.error)
            assert expected_error not in expected_errors, (
                "collect_cases merges cases into a single stubtest invocation; we already "
                "expect an error for {}".format(expected_error)
            )
            expected_errors.add(expected_error)
        output = run_stubtest(
            stub="\n\n".join(textwrap.dedent(c.stub.lstrip("\n")) for c in cases),
            runtime="\n\n".join(textwrap.dedent(c.runtime.lstrip("\n")) for c in cases),
            options=["--generate-allowlist"],
        )

        actual_errors = set(output.splitlines())
        assert actual_errors == expected_errors, output

    return test


class StubtestUnit(unittest.TestCase):
    @collect_cases
    def test_basic_good(self) -> Iterator[Case]:
        yield Case(
            stub="def f(number: int, text: str) -> None: ...",
            runtime="def f(number, text): pass",
            error=None,
        )
        yield Case(
            stub="""
            class X:
                def f(self, number: int, text: str) -> None: ...
            """,
            runtime="""
            class X:
                def f(self, number, text): pass
            """,
            error=None,
        )

    @collect_cases
    def test_types(self) -> Iterator[Case]:
        yield Case(
            stub="def mistyped_class() -> None: ...",
            runtime="class mistyped_class: pass",
            error="mistyped_class",
        )
        yield Case(
            stub="class mistyped_fn: ...", runtime="def mistyped_fn(): pass", error="mistyped_fn"
        )
        yield Case(
            stub="""
            class X:
                def mistyped_var(self) -> int: ...
            """,
            runtime="""
            class X:
                mistyped_var = 1
            """,
            error="X.mistyped_var",
        )

    @collect_cases
    def test_arg_name(self) -> Iterator[Case]:
        yield Case(
            stub="def bad(number: int, text: str) -> None: ...",
            runtime="def bad(num, text) -> None: pass",
            error="bad",
        )
        if sys.version_info >= (3, 8):
            yield Case(
                stub="def good_posonly(__number: int, text: str) -> None: ...",
                runtime="def good_posonly(num, /, text): pass",
                error=None,
            )
            yield Case(
                stub="def bad_posonly(__number: int, text: str) -> None: ...",
                runtime="def bad_posonly(flag, /, text): pass",
                error="bad_posonly",
            )
        yield Case(
            stub="""
            class BadMethod:
                def f(self, number: int, text: str) -> None: ...
            """,
            runtime="""
            class BadMethod:
                def f(self, n, text): pass
            """,
            error="BadMethod.f",
        )
        yield Case(
            stub="""
            class GoodDunder:
                def __exit__(self, t, v, tb) -> None: ...
            """,
            runtime="""
            class GoodDunder:
                def __exit__(self, exc_type, exc_val, exc_tb): pass
            """,
            error=None,
        )

    @collect_cases
    def test_arg_kind(self) -> Iterator[Case]:
        yield Case(
            stub="def runtime_kwonly(number: int, text: str) -> None: ...",
            runtime="def runtime_kwonly(number, *, text): pass",
            error="runtime_kwonly",
        )
        yield Case(
            stub="def stub_kwonly(number: int, *, text: str) -> None: ...",
            runtime="def stub_kwonly(number, text): pass",
            error="stub_kwonly",
        )
        yield Case(
            stub="def stub_posonly(__number: int, text: str) -> None: ...",
            runtime="def stub_posonly(number, text): pass",
            error="stub_posonly",
        )
        if sys.version_info >= (3, 8):
            yield Case(
                stub="def good_posonly(__number: int, text: str) -> None: ...",
                runtime="def good_posonly(number, /, text): pass",
                error=None,
            )
            yield Case(
                stub="def runtime_posonly(number: int, text: str) -> None: ...",
                runtime="def runtime_posonly(number, /, text): pass",
                error="runtime_posonly",
            )

    @collect_cases
    def test_default_value(self) -> Iterator[Case]:
        yield Case(
            stub="def f1(text: str = ...) -> None: ...",
            runtime="def f1(text = 'asdf'): pass",
            error=None,
        )
        yield Case(
            stub="def f2(text: str = ...) -> None: ...", runtime="def f2(text): pass", error="f2"
        )
        yield Case(
            stub="def f3(text: str) -> None: ...",
            runtime="def f3(text = 'asdf'): pass",
            error="f3",
        )
        yield Case(
            stub="def f4(text: str = ...) -> None: ...",
            runtime="def f4(text = None): pass",
            error="f4",
        )
        yield Case(
            stub="def f5(data: bytes = ...) -> None: ...",
            runtime="def f5(data = 'asdf'): pass",
            error="f5",
        )
        yield Case(
            stub="""
            from typing import TypeVar
            T = TypeVar("T", bound=str)
            def f6(text: T = ...) -> None: ...
            """,
            runtime="def f6(text = None): pass",
            error="f6",
        )

    @collect_cases
    def test_static_class_method(self) -> Iterator[Case]:
        yield Case(
            stub="""
            class Good:
                @classmethod
                def f(cls, number: int, text: str) -> None: ...
            """,
            runtime="""
            class Good:
                @classmethod
                def f(cls, number, text): pass
            """,
            error=None,
        )
        yield Case(
            stub="""
            class Bad1:
                def f(cls, number: int, text: str) -> None: ...
            """,
            runtime="""
            class Bad1:
                @classmethod
                def f(cls, number, text): pass
            """,
            error="Bad1.f",
        )
        yield Case(
            stub="""
            class Bad2:
                @classmethod
                def f(cls, number: int, text: str) -> None: ...
            """,
            runtime="""
            class Bad2:
                @staticmethod
                def f(self, number, text): pass
            """,
            error="Bad2.f",
        )
        yield Case(
            stub="""
            class Bad3:
                @staticmethod
                def f(cls, number: int, text: str) -> None: ...
            """,
            runtime="""
            class Bad3:
                @classmethod
                def f(self, number, text): pass
            """,
            error="Bad3.f",
        )
        yield Case(
            stub="""
            class GoodNew:
                def __new__(cls, *args, **kwargs): ...
            """,
            runtime="""
            class GoodNew:
                def __new__(cls, *args, **kwargs): pass
            """,
            error=None,
        )

    @collect_cases
    def test_arg_mismatch(self) -> Iterator[Case]:
        yield Case(
            stub="def f1(a, *, b, c) -> None: ...", runtime="def f1(a, *, b, c): pass", error=None
        )
        yield Case(
            stub="def f2(a, *, b) -> None: ...", runtime="def f2(a, *, b, c): pass", error="f2"
        )
        yield Case(
            stub="def f3(a, *, b, c) -> None: ...", runtime="def f3(a, *, b): pass", error="f3"
        )
        yield Case(
            stub="def f4(a, *, b, c) -> None: ...", runtime="def f4(a, b, *, c): pass", error="f4"
        )
        yield Case(
            stub="def f5(a, b, *, c) -> None: ...", runtime="def f5(a, *, b, c): pass", error="f5"
        )

    @collect_cases
    def test_varargs_varkwargs(self) -> Iterator[Case]:
        yield Case(
            stub="def f1(*args, **kwargs) -> None: ...",
            runtime="def f1(*args, **kwargs): pass",
            error=None,
        )
        yield Case(
            stub="def f2(*args, **kwargs) -> None: ...",
            runtime="def f2(**kwargs): pass",
            error="f2",
        )
        yield Case(
            stub="def g1(a, b, c, d) -> None: ...", runtime="def g1(a, *args): pass", error=None
        )
        yield Case(
            stub="def g2(a, b, c, d, *args) -> None: ...", runtime="def g2(a): pass", error="g2"
        )
        yield Case(
            stub="def g3(a, b, c, d, *args) -> None: ...",
            runtime="def g3(a, *args): pass",
            error=None,
        )
        yield Case(
            stub="def h1(a) -> None: ...", runtime="def h1(a, b, c, d, *args): pass", error="h1"
        )
        yield Case(
            stub="def h2(a, *args) -> None: ...", runtime="def h2(a, b, c, d): pass", error="h2"
        )
        yield Case(
            stub="def h3(a, *args) -> None: ...",
            runtime="def h3(a, b, c, d, *args): pass",
            error="h3",
        )
        yield Case(
            stub="def j1(a: int, *args) -> None: ...", runtime="def j1(a): pass", error="j1"
        )
        yield Case(
            stub="def j2(a: int) -> None: ...", runtime="def j2(a, *args): pass", error="j2"
        )
        yield Case(
            stub="def j3(a, b, c) -> None: ...", runtime="def j3(a, *args, c): pass", error="j3"
        )
        yield Case(stub="def k1(a, **kwargs) -> None: ...", runtime="def k1(a): pass", error="k1")
        yield Case(
            # In theory an error, but led to worse results in practice
            stub="def k2(a) -> None: ...",
            runtime="def k2(a, **kwargs): pass",
            error=None,
        )
        yield Case(
            stub="def k3(a, b) -> None: ...", runtime="def k3(a, **kwargs): pass", error="k3"
        )
        yield Case(
            stub="def k4(a, *, b) -> None: ...", runtime="def k4(a, **kwargs): pass", error=None
        )
        yield Case(
            stub="def k5(a, *, b) -> None: ...",
            runtime="def k5(a, *, b, c, **kwargs): pass",
            error="k5",
        )
        yield Case(
            stub="def k6(a, *, b, **kwargs) -> None: ...",
            runtime="def k6(a, *, b, c, **kwargs): pass",
            error="k6",
        )

    @collect_cases
    def test_overload(self) -> Iterator[Case]:
        yield Case(
            stub="""
            from typing import overload

            @overload
            def f1(a: int, *, c: int = ...) -> int: ...
            @overload
            def f1(a: int, b: int, c: int = ...) -> str: ...
            """,
            runtime="def f1(a, b = 0, c = 0): pass",
            error=None,
        )
        yield Case(
            stub="""
            @overload
            def f2(a: int, *, c: int = ...) -> int: ...
            @overload
            def f2(a: int, b: int, c: int = ...) -> str: ...
            """,
            runtime="def f2(a, b, c = 0): pass",
            error="f2",
        )
        yield Case(
            stub="""
            @overload
            def f3(a: int) -> int: ...
            @overload
            def f3(a: int, b: str) -> str: ...
            """,
            runtime="def f3(a, b = None): pass",
            error="f3",
        )
        yield Case(
            stub="""
            @overload
            def f4(a: int, *args, b: int, **kwargs) -> int: ...
            @overload
            def f4(a: str, *args, b: int, **kwargs) -> str: ...
            """,
            runtime="def f4(a, *args, b, **kwargs): pass",
            error=None,
        )
        if sys.version_info >= (3, 8):
            yield Case(
                stub="""
                @overload
                def f5(__a: int) -> int: ...
                @overload
                def f5(__b: str) -> str: ...
                """,
                runtime="def f5(x, /): pass",
                error=None,
            )

    @collect_cases
    def test_property(self) -> Iterator[Case]:
        yield Case(
            stub="""
            class Good:
                @property
                def f(self) -> int: ...
            """,
            runtime="""
            class Good:
                @property
                def f(self) -> int: return 1
            """,
            error=None,
        )
        yield Case(
            stub="""
            class Bad:
                @property
                def f(self) -> int: ...
            """,
            runtime="""
            class Bad:
                def f(self) -> int: return 1
            """,
            error="Bad.f",
        )
        yield Case(
            stub="""
            class GoodReadOnly:
                @property
                def f(self) -> int: ...
            """,
            runtime="""
            class GoodReadOnly:
                f = 1
            """,
            error=None,
        )
        yield Case(
            stub="""
            class BadReadOnly:
                @property
                def f(self) -> str: ...
            """,
            runtime="""
            class BadReadOnly:
                f = 1
            """,
            error="BadReadOnly.f",
        )

    @collect_cases
    def test_var(self) -> Iterator[Case]:
        yield Case(stub="x1: int", runtime="x1 = 5", error=None)
        yield Case(stub="x2: str", runtime="x2 = 5", error="x2")
        yield Case("from typing import Tuple", "", None)  # dummy case
        yield Case(
            stub="""
            x3: Tuple[int, int]
            """,
            runtime="x3 = (1, 3)",
            error=None,
        )
        yield Case(
            stub="""
            x4: Tuple[int, int]
            """,
            runtime="x4 = (1, 3, 5)",
            error="x4",
        )
        yield Case(stub="x5: int", runtime="def x5(a, b): pass", error="x5")
        yield Case(
            stub="def foo(a: int, b: int) -> None: ...\nx6 = foo",
            runtime="def foo(a, b): pass\ndef x6(c, d): pass",
            error="x6",
        )
        yield Case(
            stub="""
            class X:
                f: int
            """,
            runtime="""
            class X:
                def __init__(self):
                    self.f = "asdf"
            """,
            error=None,
        )

    @collect_cases
    def test_type_alias(self) -> Iterator[Case]:
        yield Case(
            stub="""
            class X:
                def f(self) -> None: ...
            Y = X
            """,
            runtime="""
            class X:
                def f(self) -> None: ...
            class Y: ...
            """,
            error="Y.f",
        )
        yield Case(
            stub="""
            from typing import Tuple
            A = Tuple[int, str]
            """,
            runtime="A = (int, str)",
            error="A",
        )

    @collect_cases
    def test_enum(self) -> Iterator[Case]:
        yield Case(
            stub="""
            import enum
            class X(enum.Enum):
                a: int
                b: str
                c: str
            """,
            runtime="""
            import enum
            class X(enum.Enum):
                a = 1
                b = "asdf"
                c = 2
            """,
            error="X.c",
        )

    @collect_cases
    def test_decorator(self) -> Iterator[Case]:
        yield Case(
            stub="""
            from typing import Any, Callable
            def decorator(f: Callable[[], int]) -> Callable[..., Any]: ...
            @decorator
            def f() -> Any: ...
            """,
            runtime="""
            def decorator(f): return f
            @decorator
            def f(): return 3
            """,
            error=None,
        )

    @collect_cases
    def test_missing(self) -> Iterator[Case]:
        yield Case(stub="x = 5", runtime="", error="x")
        yield Case(stub="def f(): ...", runtime="", error="f")
        yield Case(stub="class X: ...", runtime="", error="X")
        yield Case(
            stub="""
            from typing import overload
            @overload
            def h(x: int): ...
            @overload
            def h(x: str): ...
            """,
            runtime="",
            error="h",
        )
        yield Case(stub="", runtime="__all__ = []", error=None)  # dummy case
        yield Case(stub="", runtime="__all__ += ['y']\ny = 5", error="y")
        yield Case(stub="", runtime="__all__ += ['g']\ndef g(): pass", error="g")
        # Here we should only check that runtime has B, since the stub explicitly re-exports it
        yield Case(
            stub="from mystery import A, B as B, C as D  # type: ignore", runtime="", error="B"
        )

    @collect_cases
    def test_missing_no_runtime_all(self) -> Iterator[Case]:
        yield Case(stub="", runtime="import sys", error=None)
        yield Case(stub="", runtime="def g(): ...", error="g")
        yield Case(stub="", runtime="CONSTANT = 0", error="CONSTANT")

    @collect_cases
    def test_non_public_1(self) -> Iterator[Case]:
        yield Case(stub="__all__: list[str]", runtime="", error=None)  # dummy case
        yield Case(stub="_f: int", runtime="def _f(): ...", error="_f")

    @collect_cases
    def test_non_public_2(self) -> Iterator[Case]:
        yield Case(
            stub="__all__: list[str] = ['f']", runtime="__all__ = ['f']", error=None
        )
        yield Case(stub="f: int", runtime="def f(): ...", error="f")
        yield Case(stub="g: int", runtime="def g(): ...", error="g")

    @collect_cases
    def test_special_dunders(self) -> Iterator[Case]:
        yield Case(
            stub="class A:\n  def __init__(self, a: int, b: int) -> None: ...",
            runtime="class A:\n  def __init__(self, a, bx): pass",
            error="A.__init__",
        )
        yield Case(
            stub="class B:\n  def __call__(self, c: int, d: int) -> None: ...",
            runtime="class B:\n  def __call__(self, c, dx): pass",
            error="B.__call__",
        )
        if sys.version_info >= (3, 6):
            yield Case(
                stub="class C:\n  def __init_subclass__(cls, e: int, **kwargs: int) -> None: ...",
                runtime="class C:\n  def __init_subclass__(cls, e, **kwargs): pass",
                error=None,
            )
        if sys.version_info >= (3, 9):
            yield Case(
                stub="class D:\n  def __class_getitem__(cls, type: type) -> type: ...",
                runtime="class D:\n  def __class_getitem__(cls, type): ...",
                error=None,
            )

    @collect_cases
    def test_name_mangling(self) -> Iterator[Case]:
        yield Case(
            stub="""
            class X:
                def __mangle_good(self, text: str) -> None: ...
                def __mangle_bad(self, number: int) -> None: ...
            """,
            runtime="""
            class X:
                def __mangle_good(self, text): pass
                def __mangle_bad(self, text): pass
            """,
            error="X.__mangle_bad"
        )

    @collect_cases
    def test_mro(self) -> Iterator[Case]:
        yield Case(
            stub="""
            class A:
                def foo(self, x: int) -> None: ...
            class B(A):
                pass
            class C(A):
                pass
            """,
            runtime="""
            class A:
                def foo(self, x: int) -> None: ...
            class B(A):
                def foo(self, x: int) -> None: ...
            class C(A):
                def foo(self, y: int) -> None: ...
            """,
            error="C.foo"
        )
        yield Case(
            stub="""
            class X: ...
            """,
            runtime="""
            class X:
                def __init__(self, x): pass
            """,
            error="X.__init__"
        )

    @collect_cases
    def test_good_literal(self) -> Iterator[Case]:
        yield Case(
            stub=r"""
            from typing_extensions import Literal

            import enum
            class Color(enum.Enum):
                RED: int

            NUM: Literal[1]
            CHAR: Literal['a']
            FLAG: Literal[True]
            NON: Literal[None]
            BYT1: Literal[b'abc']
            BYT2: Literal[b'\x90']
            ENUM: Literal[Color.RED]
            """,
            runtime=r"""
            import enum
            class Color(enum.Enum):
                RED = 3

            NUM = 1
            CHAR = 'a'
            NON = None
            FLAG = True
            BYT1 = b"abc"
            BYT2 = b'\x90'
            ENUM = Color.RED
            """,
            error=None,
        )

    @collect_cases
    def test_bad_literal(self) -> Iterator[Case]:
        yield Case("from typing_extensions import Literal", "", None)  # dummy case
        yield Case(
            stub="INT_FLOAT_MISMATCH: Literal[1]",
            runtime="INT_FLOAT_MISMATCH = 1.0",
            error="INT_FLOAT_MISMATCH",
        )
        yield Case(
            stub="WRONG_INT: Literal[1]",
            runtime="WRONG_INT = 2",
            error="WRONG_INT",
        )
        yield Case(
            stub="WRONG_STR: Literal['a']",
            runtime="WRONG_STR = 'b'",
            error="WRONG_STR",
        )
        yield Case(
            stub="BYTES_STR_MISMATCH: Literal[b'value']",
            runtime="BYTES_STR_MISMATCH = 'value'",
            error="BYTES_STR_MISMATCH",
        )
        yield Case(
            stub="STR_BYTES_MISMATCH: Literal['value']",
            runtime="STR_BYTES_MISMATCH = b'value'",
            error="STR_BYTES_MISMATCH",
        )
        yield Case(
            stub="WRONG_BYTES: Literal[b'abc']",
            runtime="WRONG_BYTES = b'xyz'",
            error="WRONG_BYTES",
        )
        yield Case(
            stub="WRONG_BOOL_1: Literal[True]",
            runtime="WRONG_BOOL_1 = False",
            error='WRONG_BOOL_1',
        )
        yield Case(
            stub="WRONG_BOOL_2: Literal[False]",
            runtime="WRONG_BOOL_2 = True",
            error='WRONG_BOOL_2',
        )


def remove_color_code(s: str) -> str:
    return re.sub("\\x1b.*?m", "", s)  # this works!


class StubtestMiscUnit(unittest.TestCase):
    def test_output(self) -> None:
        output = run_stubtest(
            stub="def bad(number: int, text: str) -> None: ...",
            runtime="def bad(num, text): pass",
            options=[],
        )
        expected = (
            'error: {0}.bad is inconsistent, stub argument "number" differs from runtime '
            'argument "num"\nStub: at line 1\ndef (number: builtins.int, text: builtins.str)\n'
            "Runtime: at line 1 in file {0}.py\ndef (num, text)\n\n".format(TEST_MODULE_NAME)
        )
        assert remove_color_code(output) == expected

        output = run_stubtest(
            stub="def bad(number: int, text: str) -> None: ...",
            runtime="def bad(num, text): pass",
            options=["--concise"],
        )
        expected = (
            "{}.bad is inconsistent, "
            'stub argument "number" differs from runtime argument "num"\n'.format(TEST_MODULE_NAME)
        )
        assert remove_color_code(output) == expected

    def test_ignore_flags(self) -> None:
        output = run_stubtest(
            stub="", runtime="__all__ = ['f']\ndef f(): pass", options=["--ignore-missing-stub"]
        )
        assert not output

        output = run_stubtest(
            stub="", runtime="def f(): pass", options=["--ignore-missing-stub"]
        )
        assert not output

        output = run_stubtest(
            stub="def f(__a): ...", runtime="def f(a): pass", options=["--ignore-positional-only"]
        )
        assert not output

    def test_allowlist(self) -> None:
        # Can't use this as a context because Windows
        allowlist = tempfile.NamedTemporaryFile(mode="w+", delete=False)
        try:
            with allowlist:
                allowlist.write("{}.bad  # comment\n# comment".format(TEST_MODULE_NAME))

            output = run_stubtest(
                stub="def bad(number: int, text: str) -> None: ...",
                runtime="def bad(asdf, text): pass",
                options=["--allowlist", allowlist.name],
            )
            assert not output

            # test unused entry detection
            output = run_stubtest(stub="", runtime="", options=["--allowlist", allowlist.name])
            assert output == "note: unused allowlist entry {}.bad\n".format(TEST_MODULE_NAME)

            output = run_stubtest(
                stub="",
                runtime="",
                options=["--allowlist", allowlist.name, "--ignore-unused-allowlist"],
            )
            assert not output

            # test regex matching
            with open(allowlist.name, mode="w+") as f:
                f.write("{}.b.*\n".format(TEST_MODULE_NAME))
                f.write("(unused_missing)?\n")
                f.write("unused.*\n")

            output = run_stubtest(
                stub=textwrap.dedent(
                    """
                    def good() -> None: ...
                    def bad(number: int) -> None: ...
                    def also_bad(number: int) -> None: ...
                    """.lstrip("\n")
                ),
                runtime=textwrap.dedent(
                    """
                    def good(): pass
                    def bad(asdf): pass
                    def also_bad(asdf): pass
                    """.lstrip("\n")
                ),
                options=["--allowlist", allowlist.name, "--generate-allowlist"],
            )
            assert output == "note: unused allowlist entry unused.*\n{}.also_bad\n".format(
                TEST_MODULE_NAME
            )
        finally:
            os.unlink(allowlist.name)

    def test_mypy_build(self) -> None:
        output = run_stubtest(stub="+", runtime="", options=[])
        assert remove_color_code(output) == (
            "error: not checking stubs due to failed mypy compile:\n{}.pyi:1: "
            "error: invalid syntax\n".format(TEST_MODULE_NAME)
        )

        output = run_stubtest(stub="def f(): ...\ndef f(): ...", runtime="", options=[])
        assert remove_color_code(output) == (
            'error: not checking stubs due to mypy build errors:\n{}.pyi:2: '
            'error: Name "f" already defined on line 1\n'.format(TEST_MODULE_NAME)
        )

    def test_missing_stubs(self) -> None:
        output = io.StringIO()
        with contextlib.redirect_stdout(output):
            test_stubs(parse_options(["not_a_module"]))
        assert "error: not_a_module failed to find stubs" in remove_color_code(output.getvalue())

    def test_get_typeshed_stdlib_modules(self) -> None:
        stdlib = mypy.stubtest.get_typeshed_stdlib_modules(None)
        assert "builtins" in stdlib
        assert "os" in stdlib
        assert "os.path" in stdlib
        assert "asyncio" in stdlib
        assert ("dataclasses" in stdlib) == (sys.version_info >= (3, 7))

    def test_signature(self) -> None:
        def f(a: int, b: int, *, c: int, d: int = 0, **kwargs: Any) -> None:
            pass

        assert (
            str(mypy.stubtest.Signature.from_inspect_signature(inspect.signature(f)))
            == "def (a, b, *, c, d = ..., **kwargs)"
        )

    def test_config_file(self) -> None:
        runtime = "temp = 5\n"
        stub = "from decimal import Decimal\ntemp: Decimal\n"
        config_file = (
            "[mypy]\n"
            "plugins={}/test-data/unit/plugins/decimal_to_int.py\n".format(root_dir)
        )
        output = run_stubtest(stub=stub, runtime=runtime, options=[])
        assert output == (
            "error: test_module.temp variable differs from runtime type Literal[5]\n"
            "Stub: at line 2\ndecimal.Decimal\nRuntime:\n5\n\n"
        )
        output = run_stubtest(stub=stub, runtime=runtime, options=[], config_file=config_file)
        assert output == ""
