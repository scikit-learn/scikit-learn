"""Test that C functions used in primitives are declared in a header such as CPy.h."""

import glob
import os
import re
import unittest

from mypyc.primitives import registry
from mypyc.primitives.registry import CFunctionDescription


class TestHeaderInclusion(unittest.TestCase):
    def test_primitives_included_in_header(self) -> None:
        base_dir = os.path.join(os.path.dirname(__file__), '..', 'lib-rt')
        with open(os.path.join(base_dir, 'CPy.h')) as f:
            header = f.read()
        with open(os.path.join(base_dir, 'pythonsupport.h')) as f:
            header += f.read()

        def check_name(name: str) -> None:
            if name.startswith('CPy'):
                assert re.search(r'\b{}\b'.format(name), header), (
                    '"{}" is used in mypyc.primitives but not declared in CPy.h'.format(name))

        for values in [registry.method_call_ops.values(),
                       registry.function_ops.values(),
                       registry.binary_ops.values(),
                       registry.unary_ops.values()]:
            for ops in values:
                if isinstance(ops, CFunctionDescription):
                    ops = [ops]
                for op in ops:
                    check_name(op.c_function_name)

        primitives_path = os.path.join(os.path.dirname(__file__), '..', 'primitives')
        for fnam in glob.glob('{}/*.py'.format(primitives_path)):
            with open(fnam) as f:
                content = f.read()
            for name in re.findall(r'c_function_name=["\'](CPy[A-Z_a-z0-9]+)', content):
                check_name(name)
