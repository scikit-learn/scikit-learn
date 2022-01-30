from typing import Dict, List, Sequence

cmp_op: Sequence[str]
hasconst: List[int]
hasname: List[int]
hasjrel: List[int]
hasjabs: List[int]
haslocal: List[int]
hascompare: List[int]
hasfree: List[int]
opname: List[str]

opmap: Dict[str, int]
HAVE_ARGUMENT: int
EXTENDED_ARG: int
