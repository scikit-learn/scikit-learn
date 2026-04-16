"""The analyses submodule contains all the analyses passes offered in Pythran.

This file is just for convenience and turns the import from

import analyses.foo.Foo

into

import analyses.Foo
"""

from .aliases import Aliases, StrictAliases, InterproceduralAliases
from .ancestors import Ancestors, AncestorsWithBody
from .argument_effects import ArgumentEffects
from .argument_read_once import ArgumentReadOnce
from .ast_matcher import ASTMatcher, AST_any, AST_or, Placeholder, Check
from .cfg import CFG
from .constant_expressions import ConstantExpressions
from .dependencies import Dependencies
from .extended_syntax_check import ExtendedSyntaxCheck
from .fixed_size_list import FixedSizeList
from .global_declarations import GlobalDeclarations, NonlocalDeclarations
from .global_effects import GlobalEffects
from .globals_analysis import Globals
from .has_return import HasReturn, HasBreak, HasContinue
from .identifiers import Identifiers
from .immediates import Immediates
from .imported_ids import ImportedIds
from .inlinable import Inlinable
from .is_assigned import IsAssigned
from .lazyness_analysis import LazynessAnalysis
from .literals import Literals
from .local_declarations import LocalNodeDeclarations, LocalNameDeclarations
from .locals_analysis import Locals
from .node_count import NodeCount
from .optimizable_comprehension import OptimizableComprehension
from .ordered_global_declarations import OrderedGlobalDeclarations
from .parallel_maps import ParallelMaps
from .potential_iterator import PotentialIterator
from .pure_expressions import PureExpressions
from .pure_functions import PureFunctions
from .range_values import RangeValues
from .scope import Scope
from .static_expressions import StaticExpressions, HasStaticExpression
from .use_def_chain import DefUseChains, UseDefChains, ExtendedDefUseChains
from .use_omp import UseOMP
from .yield_points import YieldPoints
