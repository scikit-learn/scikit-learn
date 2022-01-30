import unittest

from mypyc.codegen.emitclass import slot_key


class TestEmitClass(unittest.TestCase):
    def test_slot_key(self) -> None:
        attrs = ['__add__', '__radd__', '__rshift__', '__rrshift__', '__setitem__', '__delitem__']
        s = sorted(attrs, key=lambda x: slot_key(x))
        # __delitem__ and reverse methods should come last.
        assert s == [
            '__add__', '__rshift__', '__setitem__', '__delitem__', '__radd__', '__rrshift__']
