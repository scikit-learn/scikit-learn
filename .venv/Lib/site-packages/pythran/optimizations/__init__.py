"""
The optimisations submodule contains all the optimisations offered in Pythran.

This file is just for convenience and turns the import from

import optimisations.xxxxx.xxxxx

into

import optimisations.xxxxx
"""

from .constant_folding import ConstantFolding, PartialConstantFolding
from .copyto import CopyTo
from .dead_code_elimination import DeadCodeElimination
from .fast_gexpr import FastGExpr
from .forward_substitution import ForwardSubstitution, PreInliningForwardSubstitution
from .iter_transformation import IterTransformation
from .comprehension_patterns import ComprehensionPatterns
from .list_comp_to_genexp import ListCompToGenexp
from .loop_full_unrolling import LoopFullUnrolling
from .modindex import ModIndex
from .pattern_transform import PatternTransform
from .range_loop_unfolding import RangeLoopUnfolding
from .range_based_simplify import RangeBasedSimplify
from .square import Square
from .inlining import Inlining
from .inline_builtins import InlineBuiltins
from .list_to_tuple import ListToTuple
from .tuple_to_shape import TupleToShape
from .remove_dead_functions import RemoveDeadFunctions
from .simplify_except import SimplifyExcept
