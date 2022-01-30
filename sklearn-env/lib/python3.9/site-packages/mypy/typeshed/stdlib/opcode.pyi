import sys
from typing import Sequence

cmp_op: Sequence[str]
hasconst: list[int]
hasname: list[int]
hasjrel: list[int]
hasjabs: list[int]
haslocal: list[int]
hascompare: list[int]
hasfree: list[int]
opname: list[str]

opmap: dict[str, int]
HAVE_ARGUMENT: int
EXTENDED_ARG: int

if sys.version_info >= (3, 8):
    def stack_effect(__opcode: int, __oparg: int | None = ..., *, jump: bool | None = ...) -> int: ...

else:
    def stack_effect(__opcode: int, __oparg: int | None = ...) -> int: ...

hasnargs: list[int]
