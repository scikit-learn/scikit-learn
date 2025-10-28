# SPDX-License-Identifier: Apache-2.0
# Copyright 2019 The Meson development team

from __future__ import annotations

from .. import mesonlib
from .. import mlog
from .common import cmake_is_debug
import typing as T

if T.TYPE_CHECKING:
    from .traceparser import CMakeTraceParser, CMakeTarget

def parse_generator_expressions(
            raw: str,
            trace: 'CMakeTraceParser',
            *,
            context_tgt: T.Optional['CMakeTarget'] = None,
        ) -> str:
    '''Parse CMake generator expressions

    Most generator expressions are simply ignored for
    simplicity, however some are required for some common
    use cases.
    '''

    # Early abort if no generator expression present
    if '$<' not in raw:
        return raw

    out = ''
    i = 0

    def equal(arg: str) -> str:
        col_pos = arg.find(',')
        if col_pos < 0:
            return '0'
        else:
            return '1' if arg[:col_pos] == arg[col_pos + 1:] else '0'

    def vers_comp(op: str, arg: str) -> str:
        col_pos = arg.find(',')
        if col_pos < 0:
            return '0'
        else:
            return '1' if mesonlib.version_compare(arg[:col_pos], '{}{}'.format(op, arg[col_pos + 1:])) else '0'

    def target_property(arg: str) -> str:
        # We can't really support this since we don't have any context
        if ',' not in arg:
            if context_tgt is None:
                return ''
            return ';'.join(context_tgt.properties.get(arg, []))

        args = arg.split(',')
        props = trace.targets[args[0]].properties.get(args[1], []) if args[0] in trace.targets else []
        return ';'.join(props)

    def target_file(arg: str) -> str:
        if arg not in trace.targets:
            mlog.warning(f"Unable to evaluate the cmake variable '$<TARGET_FILE:{arg}>'.")
            return ''
        tgt = trace.targets[arg]

        cfgs = []
        cfg = ''

        if 'IMPORTED_CONFIGURATIONS' in tgt.properties:
            cfgs = [x for x in tgt.properties['IMPORTED_CONFIGURATIONS'] if x]
            cfg = cfgs[0]

        if cmake_is_debug(trace.env):
            if 'DEBUG' in cfgs:
                cfg = 'DEBUG'
            elif 'RELEASE' in cfgs:
                cfg = 'RELEASE'
        else:
            if 'RELEASE' in cfgs:
                cfg = 'RELEASE'

        if f'IMPORTED_IMPLIB_{cfg}' in tgt.properties:
            return ';'.join([x for x in tgt.properties[f'IMPORTED_IMPLIB_{cfg}'] if x])
        elif 'IMPORTED_IMPLIB' in tgt.properties:
            return ';'.join([x for x in tgt.properties['IMPORTED_IMPLIB'] if x])
        elif f'IMPORTED_LOCATION_{cfg}' in tgt.properties:
            return ';'.join([x for x in tgt.properties[f'IMPORTED_LOCATION_{cfg}'] if x])
        elif 'IMPORTED_LOCATION' in tgt.properties:
            return ';'.join([x for x in tgt.properties['IMPORTED_LOCATION'] if x])
        return ''

    supported: T.Dict[str, T.Callable[[str], str]] = {
        # Boolean functions
        'BOOL': lambda x: '0' if x.upper() in {'', '0', 'FALSE', 'OFF', 'N', 'NO', 'IGNORE', 'NOTFOUND'} or x.endswith('-NOTFOUND') else '1',
        'AND': lambda x: '1' if all(y == '1' for y in x.split(',')) else '0',
        'OR': lambda x: '1' if any(y == '1' for y in x.split(',')) else '0',
        'NOT': lambda x: '0' if x == '1' else '1',

        'IF': lambda x: x.split(',')[1] if x.split(',')[0] == '1' else x.split(',')[2],

        '0': lambda x: '',
        '1': lambda x: x,

        # String operations
        'STREQUAL': equal,
        'EQUAL': equal,
        'VERSION_LESS': lambda x: vers_comp('<', x),
        'VERSION_GREATER': lambda x: vers_comp('>', x),
        'VERSION_EQUAL': lambda x: vers_comp('=', x),
        'VERSION_LESS_EQUAL': lambda x: vers_comp('<=', x),
        'VERSION_GREATER_EQUAL': lambda x: vers_comp('>=', x),

        # String modification
        'LOWER_CASE': lambda x: x.lower(),
        'UPPER_CASE': lambda x: x.upper(),

        # Always assume the BUILD_INTERFACE is valid.
        # INSTALL_INTERFACE is always invalid for subprojects and
        # it should also never appear in CMake config files, used
        # for dependencies
        'INSTALL_INTERFACE': lambda x: '',
        'BUILD_INTERFACE': lambda x: x,

        # Constants
        'ANGLE-R': lambda x: '>',
        'COMMA': lambda x: ',',
        'SEMICOLON': lambda x: ';',

        # Target related expressions
        'TARGET_EXISTS': lambda x: '1' if x in trace.targets else '0',
        'TARGET_NAME_IF_EXISTS': lambda x: x if x in trace.targets else '',
        'TARGET_PROPERTY': target_property,
        'TARGET_FILE': target_file,
    }

    # Recursively evaluate generator expressions
    def eval_generator_expressions() -> str:
        nonlocal i
        i += 2

        func = ''
        args = ''
        res = ''
        exp = ''

        # Determine the body of the expression
        while i < len(raw):
            if raw[i] == '>':
                # End of the generator expression
                break
            elif i < len(raw) - 1 and raw[i] == '$' and raw[i + 1] == '<':
                # Nested generator expression
                exp += eval_generator_expressions()
            else:
                # Generator expression body
                exp += raw[i]

            i += 1

        # Split the expression into a function and arguments part
        col_pos = exp.find(':')
        if col_pos < 0:
            func = exp
        else:
            func = exp[:col_pos]
            args = exp[col_pos + 1:]

        func = func.strip()
        args = args.strip()

        # Evaluate the function
        if func in supported:
            res = supported[func](args)

        return res

    while i < len(raw):
        if i < len(raw) - 1 and raw[i] == '$' and raw[i + 1] == '<':
            # Generator expression detected --> try resolving it
            out += eval_generator_expressions()
        else:
            # Normal string, leave unchanged
            out += raw[i]

        i += 1

    return out
