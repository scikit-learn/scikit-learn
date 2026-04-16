# SPDX-License-Identifier: Apache-2.0

from enum import Enum
import typing as T

class MesonOperator(Enum):
    # Arithmetic
    PLUS = '+'
    MINUS = '-'
    TIMES = '*'
    DIV = '/'
    MOD = '%'

    UMINUS = 'uminus'

    # Logic
    NOT = 'not'

    # Should return the boolsche interpretation of the value (`'' == false` for instance)
    BOOL = 'bool()'

    # Comparison
    EQUALS = '=='
    NOT_EQUALS = '!='
    GREATER = '>'
    LESS = '<'
    GREATER_EQUALS = '>='
    LESS_EQUALS = '<='

    # Container
    IN = 'in'
    NOT_IN = 'not in'
    INDEX = '[]'

# Accessing this directly is about 9x faster than calling MesonOperator(s),
# and about 3 times faster than a staticmethod
MAPPING: T.Mapping[str, MesonOperator] = {x.value: x for x in MesonOperator}
