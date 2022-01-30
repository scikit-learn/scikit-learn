import unittest
from typing import Dict

from mypyc.codegen.emit import Emitter, EmitterContext
from mypyc.ir.ops import BasicBlock, Value, Register
from mypyc.ir.rtypes import int_rprimitive
from mypyc.namegen import NameGenerator


class TestEmitter(unittest.TestCase):
    def setUp(self) -> None:
        self.n = Register(int_rprimitive, 'n')
        self.context = EmitterContext(NameGenerator([['mod']]))

    def test_label(self) -> None:
        emitter = Emitter(self.context, {})
        assert emitter.label(BasicBlock(4)) == 'CPyL4'

    def test_reg(self) -> None:
        names: Dict[Value, str] = {self.n: "n"}
        emitter = Emitter(self.context, names)
        assert emitter.reg(self.n) == 'cpy_r_n'

    def test_emit_line(self) -> None:
        emitter = Emitter(self.context, {})
        emitter.emit_line('line;')
        emitter.emit_line('a {')
        emitter.emit_line('f();')
        emitter.emit_line('}')
        assert emitter.fragments == ['line;\n',
                                     'a {\n',
                                     '    f();\n',
                                     '}\n']
