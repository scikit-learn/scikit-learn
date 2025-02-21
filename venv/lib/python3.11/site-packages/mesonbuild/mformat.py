# SPDX-License-Identifier: Apache-2.0
# Copyright 2024 The Meson development team

from __future__ import annotations

import re
import typing as T
from configparser import ConfigParser, MissingSectionHeaderError, ParsingError
from copy import deepcopy
from dataclasses import dataclass, field, fields, asdict
from pathlib import Path
import sys

from . import mparser
from .mesonlib import MesonException
from .ast.postprocess import AstConditionLevel
from .ast.printer import RawPrinter
from .ast.visitor import FullAstVisitor
from .environment import build_filename

if T.TYPE_CHECKING:
    import argparse
    from typing_extensions import Literal


class DefaultConfigParser(ConfigParser):

    def __init__(self, delimiters: T.Tuple[str, ...] = ('=', ':')):
        super().__init__(delimiters=delimiters, interpolation=None)

    def read_default(self, filename: Path) -> None:
        if not filename.exists():
            raise MesonException(f'Configuration file {filename} not found')
        try:
            super().read(filename, encoding='utf-8')
        except MissingSectionHeaderError:
            self.read_string(f'[{self.default_section}]\n' + filename.read_text(encoding='utf-8'))

    def getstr(self, section: str, key: str, fallback: T.Optional[str] = None) -> T.Optional[str]:
        value: T.Optional[str] = self.get(section, key, fallback=fallback)
        if value:
            value = value.strip('"').strip("'")
        return value


def match_path(filename: str, pattern: str) -> bool:
    '''recursive glob match for editorconfig sections'''
    index = 0
    num_ranges: T.List[T.Tuple[int, int]] = []

    def curl_replace(m: re.Match) -> str:
        nonlocal index

        if '\\.\\.' in m[1]:
            index += 1
            low, high = m[1].split('\\.\\.')
            num_ranges.append((int(low), int(high)))
            return f'(?P<num{index}>-?[0-9]+)'
        else:
            return T.cast(str, m[1].replace(',', '|'))

    pattern_re = pattern.replace('.', '\\.')
    pattern_re = re.sub(r'(?<!\\)\?', '.', pattern_re)  # ? -> .
    pattern_re = re.sub(r'(?<![\\\*])\*(?!\*)', '([^/]*)', pattern_re)  # * -> ([^/]*)
    pattern_re = re.sub(r'(?<!\\)\*\*', '(.*)', pattern_re)  # ** -> (.*)
    pattern_re = re.sub(r'(?<!\\)\[!(.*?[^\\])\]', r'([^\1])', pattern_re)  # [!name] -> [^name]
    pattern_re = re.sub(r'(?<!\\)\{(.*?[^\\])}', curl_replace, pattern_re)  # {}
    if pattern.startswith('/'):
        pattern_re = '^' + pattern_re
    pattern_re += '$'

    m = re.search(pattern_re, filename)
    if m is None:
        return False

    for i in range(index):
        try:
            val = int(m[f'num{i+1}'])
            if not num_ranges[i][0] <= val <= num_ranges[i][1]:
                return False
        except ValueError:
            return False

    return True


@dataclass
class EditorConfig:

    indent_style: T.Optional[Literal['space', 'tab']] = field(default=None, metadata={'getter': DefaultConfigParser.get})
    indent_size: T.Optional[int] = field(default=None, metadata={'getter': DefaultConfigParser.getint})
    tab_width: T.Optional[int] = field(default=None, metadata={'getter': DefaultConfigParser.getint})
    end_of_line: T.Optional[Literal['lf', 'cr', 'crlf']] = field(default=None, metadata={'getter': DefaultConfigParser.get})
    charset: T.Optional[Literal['latin1', 'utf-8', 'utf-8-bom', 'utf-16be', 'utf-16le']] = field(default=None, metadata={'getter': DefaultConfigParser.get})
    trim_trailing_whitespace: T.Optional[bool] = field(default=None, metadata={'getter': DefaultConfigParser.getboolean})
    insert_final_newline: T.Optional[bool] = field(default=None, metadata={'getter': DefaultConfigParser.getboolean})
    max_line_length: T.Optional[T.Union[Literal['off'], int]] = field(default=None, metadata={'getter': DefaultConfigParser.get})


@dataclass
class FormatterConfig:

    # Config keys compatible with muon
    max_line_length: T.Optional[int] = field(
        default=None,
        metadata={'getter': DefaultConfigParser.getint,
                  'default': 80,
                  })
    indent_by: T.Optional[str] = field(
        default=None,
        metadata={'getter': DefaultConfigParser.getstr,
                  'default': '    ',
                  })
    space_array: T.Optional[bool] = field(
        default=None,
        metadata={'getter': DefaultConfigParser.getboolean,
                  'default': False,
                  })
    kwargs_force_multiline: T.Optional[bool] = field(
        default=None,  # kwa_ml
        metadata={'getter': DefaultConfigParser.getboolean,
                  'default': False,
                  })
    wide_colon: T.Optional[bool] = field(
        default=None,
        metadata={'getter': DefaultConfigParser.getboolean,
                  'default': False,
                  })
    no_single_comma_function: T.Optional[bool] = field(
        default=None,
        metadata={'getter': DefaultConfigParser.getboolean,
                  'default': False,
                  })

    # Additional config keys
    end_of_line: T.Optional[Literal['cr', 'lf', 'crlf', 'native']] = field(
        default=None,
        metadata={'getter': DefaultConfigParser.getstr,
                  'default': 'native',
                  })
    indent_before_comments: T.Optional[str] = field(
        default=None,
        metadata={'getter': DefaultConfigParser.getstr,
                  'default': '  ',
                  })
    simplify_string_literals: T.Optional[bool] = field(
        default=None,
        metadata={'getter': DefaultConfigParser.getboolean,
                  'default': True,
                  })
    insert_final_newline: T.Optional[bool] = field(
        default=None,
        metadata={'getter': DefaultConfigParser.getboolean,
                  'default': True,
                  })
    tab_width: T.Optional[int] = field(
        default=None,
        metadata={'getter': DefaultConfigParser.getint,
                  'default': 4,
                  }
    )
    sort_files: T.Optional[bool] = field(
        default=None,
        metadata={'getter': DefaultConfigParser.getboolean,
                  'default': True,
                  })
    group_arg_value: T.Optional[bool] = field(
        default=None,
        metadata={'getter': DefaultConfigParser.getboolean,
                  'default': False,
                  })
    use_editor_config: T.Optional[bool] = field(
        default=None,
        metadata={'getter': DefaultConfigParser.getboolean,
                  'default': False,
                  })

    @classmethod
    def default(cls) -> FormatterConfig:
        defaults = {f.name: f.metadata['default'] for f in fields(cls)}
        return cls(**defaults)

    def update(self, config: FormatterConfig) -> FormatterConfig:
        """Returns copy of self updated with other config"""
        new_config = deepcopy(self)
        for key, value in asdict(config).items():
            if value is not None:
                setattr(new_config, key, value)
        return new_config

    def with_editorconfig(self, editorconfig: EditorConfig) -> FormatterConfig:
        """Returns copy of self updated with editorconfig"""
        config = deepcopy(self)

        if editorconfig.indent_style == 'space':
            indent_size = editorconfig.indent_size or 4
            config.indent_by = indent_size * ' '
        elif editorconfig.indent_style == 'tab':
            config.indent_by = '\t'
        elif editorconfig.indent_size:
            config.indent_by = editorconfig.indent_size * ' '

        if editorconfig.max_line_length == 'off':
            config.max_line_length = 0
        elif editorconfig.max_line_length:
            config.max_line_length = int(editorconfig.max_line_length)

        if editorconfig.end_of_line:
            config.end_of_line = editorconfig.end_of_line
        if editorconfig.insert_final_newline:
            config.insert_final_newline = editorconfig.insert_final_newline
        if editorconfig.tab_width:
            config.tab_width = editorconfig.tab_width

        return config

    @property
    def newline(self) -> T.Optional[str]:
        if self.end_of_line == 'crlf':
            return '\r\n'
        if self.end_of_line == 'lf':
            return '\n'
        if self.end_of_line == 'cr':
            return '\r'
        return None


class MultilineArgumentDetector(FullAstVisitor):

    def __init__(self, config: FormatterConfig):
        self.config = config
        self.is_multiline = False

    def enter_node(self, node: mparser.BaseNode) -> None:
        if node.whitespaces and '#' in node.whitespaces.value:
            self.is_multiline = True

        elif isinstance(node, mparser.StringNode) and node.is_multiline:
            self.is_multiline = True

    def visit_ArgumentNode(self, node: mparser.ArgumentNode) -> None:
        if node.is_multiline:
            self.is_multiline = True

        if self.is_multiline:
            return

        if self.config.kwargs_force_multiline and node.kwargs:
            self.is_multiline = True

        super().visit_ArgumentNode(node)


class TrimWhitespaces(FullAstVisitor):

    def __init__(self, config: FormatterConfig):
        self.config = config

        self.in_block_comments = False
        self.in_arguments = 0
        self.indent_comments = ''

    def visit_default_func(self, node: mparser.BaseNode) -> None:
        self.enter_node(node)
        node.whitespaces.accept(self)

    def enter_node(self, node: mparser.BaseNode) -> None:
        if isinstance(node, mparser.WhitespaceNode):
            return
        if not node.whitespaces:
            # Ensure every node has a whitespace node
            node.whitespaces = mparser.WhitespaceNode(mparser.Token('whitespace', node.filename, 0, 0, 0, (0, 0), ''))
            node.whitespaces.condition_level = node.condition_level

    def exit_node(self, node: mparser.BaseNode) -> None:
        pass

    def move_whitespaces(self, from_node: mparser.BaseNode, to_node: mparser.BaseNode) -> None:
        to_node.whitespaces.value = from_node.whitespaces.value + to_node.whitespaces.value
        to_node.whitespaces.is_continuation = from_node.whitespaces.is_continuation
        from_node.whitespaces = None
        to_node.whitespaces.accept(self)

    def add_space_after(self, node: mparser.BaseNode) -> None:
        if not node.whitespaces.value:
            node.whitespaces.value = ' '

    def add_nl_after(self, node: mparser.BaseNode, force: bool = False) -> None:
        if not node.whitespaces.value:
            node.whitespaces.value = '\n'
        elif force and not node.whitespaces.value.endswith('\n'):
            node.whitespaces.value += '\n'

    def dedent(self, value: str) -> str:
        if value.endswith(self.config.indent_by):
            value = value[:-len(self.config.indent_by)]
        return value

    def sort_arguments(self, node: mparser.ArgumentNode) -> None:
        # TODO: natsort
        def sort_key(arg: mparser.BaseNode) -> str:
            if isinstance(arg, mparser.StringNode):
                return arg.raw_value
            return getattr(node, 'value', '')

        node.arguments.sort(key=sort_key)

    def visit_EmptyNode(self, node: mparser.EmptyNode) -> None:
        self.enter_node(node)
        self.in_block_comments = True
        node.whitespaces.accept(self)
        self.in_block_comments = False

    def visit_WhitespaceNode(self, node: mparser.WhitespaceNode) -> None:
        lines = node.value.splitlines(keepends=True)
        node.value = ''
        in_block_comments = self.in_block_comments
        with_comments = ['#' in line for line in lines] + [False]
        for i, line in enumerate(lines):
            has_nl = line.endswith('\n')
            line = line.strip()
            if line.startswith('\\'):
                node.value += ' '  # add space before \
                node.is_continuation = True
            elif line.startswith('#'):
                if not in_block_comments:
                    node.value += self.config.indent_before_comments
                else:
                    node.value += self.indent_comments
            node.value += line
            if has_nl and (line or with_comments[i+1] or not self.in_arguments):
                node.value += '\n'
            in_block_comments = True
        if node.value.endswith('\n'):
            node.value += self.indent_comments
            if node.is_continuation:
                node.value += self.config.indent_by

    def visit_SymbolNode(self, node: mparser.SymbolNode) -> None:
        super().visit_SymbolNode(node)
        if node.value in "([{" and node.whitespaces.value == '\n':
            node.whitespaces.value = ''

    def visit_StringNode(self, node: mparser.StringNode) -> None:
        self.enter_node(node)

        if self.config.simplify_string_literals:
            if node.is_multiline and not any(x in node.value for x in ['\n', "'"]):
                node.is_multiline = False
                node.value = node.escape()

            if node.is_fstring and '@' not in node.value:
                node.is_fstring = False

        self.exit_node(node)

    def visit_UnaryOperatorNode(self, node: mparser.UnaryOperatorNode) -> None:
        super().visit_UnaryOperatorNode(node)
        self.move_whitespaces(node.value, node)

    def visit_NotNode(self, node: mparser.NotNode) -> None:
        super().visit_UnaryOperatorNode(node)
        if not node.operator.whitespaces.value:
            node.operator.whitespaces.value = ' '
        self.move_whitespaces(node.value, node)

    def visit_BinaryOperatorNode(self, node: mparser.BinaryOperatorNode) -> None:
        super().visit_BinaryOperatorNode(node)
        self.add_space_after(node.left)
        self.add_space_after(node.operator)
        self.move_whitespaces(node.right, node)

    def visit_ArrayNode(self, node: mparser.ArrayNode) -> None:
        super().visit_ArrayNode(node)
        self.move_whitespaces(node.rbracket, node)

        if node.lbracket.whitespaces.value:
            node.args.is_multiline = True
        if node.args.arguments and not node.args.is_multiline and self.config.space_array:
            self.add_space_after(node.lbracket)
            self.add_space_after(node.args)
        if not node.args.arguments:
            self.move_whitespaces(node.lbracket, node.args)

    def visit_DictNode(self, node: mparser.DictNode) -> None:
        super().visit_DictNode(node)
        self.move_whitespaces(node.rcurl, node)

        if node.lcurl.whitespaces.value:
            node.args.is_multiline = True

    def visit_CodeBlockNode(self, node: mparser.CodeBlockNode) -> None:
        self.enter_node(node)
        if node.pre_whitespaces:
            self.in_block_comments = True
            node.pre_whitespaces.accept(self)
            self.in_block_comments = False
        else:
            node.pre_whitespaces = mparser.WhitespaceNode(mparser.Token('whitespace', node.filename, 0, 0, 0, (0, 0), ''))
        node.pre_whitespaces.block_indent = True

        for i in node.lines:
            i.accept(self)
        self.exit_node(node)

        if node.lines:
            self.move_whitespaces(node.lines[-1], node)
        else:
            node.whitespaces.value = node.pre_whitespaces.value + node.whitespaces.value
            node.pre_whitespaces.value = ''
            self.in_block_comments = True
            node.whitespaces.accept(self)
            self.in_block_comments = False

        if node.condition_level == 0 and self.config.insert_final_newline:
            self.add_nl_after(node, force=True)

        indent = node.condition_level * self.config.indent_by
        if indent and node.lines:
            node.pre_whitespaces.value += indent
        for line in node.lines[:-1]:
            line.whitespaces.value += indent

    def visit_IndexNode(self, node: mparser.IndexNode) -> None:
        super().visit_IndexNode(node)
        self.move_whitespaces(node.rbracket, node)

    def visit_MethodNode(self, node: mparser.MethodNode) -> None:
        super().visit_MethodNode(node)
        self.move_whitespaces(node.rpar, node)

        if node.lpar.whitespaces.value:
            node.args.is_multiline = True

    def visit_FunctionNode(self, node: mparser.FunctionNode) -> None:
        if node.func_name.value == 'files':
            if self.config.sort_files:
                self.sort_arguments(node.args)

            if len(node.args.arguments) == 1 and not node.args.kwargs:
                arg = node.args.arguments[0]
                if isinstance(arg, mparser.ArrayNode):
                    if not arg.lbracket.whitespaces or not arg.lbracket.whitespaces.value.strip():
                        # files([...]) -> files(...)
                        node.args = arg.args

        super().visit_FunctionNode(node)
        self.move_whitespaces(node.rpar, node)

        if node.lpar.whitespaces.value:
            node.args.is_multiline = True

    def visit_AssignmentNode(self, node: mparser.AssignmentNode) -> None:
        super().visit_AssignmentNode(node)
        self.add_space_after(node.var_name)
        self.add_space_after(node.operator)
        self.move_whitespaces(node.value, node)

    def visit_ForeachClauseNode(self, node: mparser.ForeachClauseNode) -> None:
        super().visit_ForeachClauseNode(node)
        self.add_space_after(node.foreach_)
        self.add_space_after(node.varnames[-1])
        for comma in node.commas:
            self.add_space_after(comma)
        self.add_space_after(node.colon)

        node.block.whitespaces.value += node.condition_level * self.config.indent_by
        node.block.whitespaces.block_indent = True

        self.move_whitespaces(node.endforeach, node)

    def visit_IfClauseNode(self, node: mparser.IfClauseNode) -> None:
        super().visit_IfClauseNode(node)
        self.move_whitespaces(node.endif, node)

        for if_node in node.ifs:
            if_node.whitespaces.value += node.condition_level * self.config.indent_by
        if isinstance(node.elseblock, mparser.ElseNode):
            node.elseblock.whitespaces.value += node.condition_level * self.config.indent_by

    def visit_IfNode(self, node: mparser.IfNode) -> None:
        super().visit_IfNode(node)
        self.add_space_after(node.if_)
        self.in_block_comments = True
        self.move_whitespaces(node.block, node)
        self.in_block_comments = False
        node.whitespaces.condition_level = node.condition_level + 1
        node.whitespaces.block_indent = True

    def visit_ElseNode(self, node: mparser.ElseNode) -> None:
        super().visit_ElseNode(node)
        self.in_block_comments = True
        self.move_whitespaces(node.block, node)
        self.in_block_comments = False
        node.whitespaces.condition_level = node.condition_level + 1
        node.whitespaces.block_indent = True

    def visit_TernaryNode(self, node: mparser.TernaryNode) -> None:
        super().visit_TernaryNode(node)
        self.add_space_after(node.condition)
        self.add_space_after(node.questionmark)
        self.add_space_after(node.trueblock)
        self.add_space_after(node.colon)
        self.move_whitespaces(node.falseblock, node)

    def visit_ArgumentNode(self, node: mparser.ArgumentNode) -> None:
        if not node.is_multiline:
            ml_detector = MultilineArgumentDetector(self.config)
            node.accept(ml_detector)
            if ml_detector.is_multiline:
                node.is_multiline = True

        self.in_arguments += 1
        super().visit_ArgumentNode(node)
        self.in_arguments -= 1

        if not node.arguments and not node.kwargs:
            node.whitespaces.accept(self)
            return

        last_node: mparser.BaseNode
        has_trailing_comma = len(node.commas) == len(node.arguments) + len(node.kwargs)
        if has_trailing_comma:
            last_node = node.commas[-1]
        elif node.kwargs:
            for last_node in node.kwargs.values():
                pass
        else:
            last_node = node.arguments[-1]

        self.move_whitespaces(last_node, node)

        if not node.is_multiline and '#' not in node.whitespaces.value:
            node.whitespaces.value = ''

    def visit_ParenthesizedNode(self, node: mparser.ParenthesizedNode) -> None:
        self.enter_node(node)

        is_multiline = node.lpar.whitespaces and '#' in node.lpar.whitespaces.value
        if is_multiline:
            self.indent_comments += self.config.indent_by

        node.lpar.accept(self)
        node.inner.accept(self)

        if is_multiline:
            node.inner.whitespaces.value = self.dedent(node.inner.whitespaces.value)
            self.indent_comments = self.dedent(self.indent_comments)
            self.add_nl_after(node.inner)

        node.rpar.accept(self)
        self.move_whitespaces(node.rpar, node)


class ArgumentFormatter(FullAstVisitor):

    def __init__(self, config: FormatterConfig):
        self.config = config
        self.level = 0
        self.indent_after = False
        self.is_function_arguments = False

    def add_space_after(self, node: mparser.BaseNode) -> None:
        if not node.whitespaces.value:
            node.whitespaces.value = ' '

    def add_nl_after(self, node: mparser.BaseNode, indent: int) -> None:
        if not node.whitespaces.value or node.whitespaces.value == ' ':
            node.whitespaces.value = '\n'
        indent_by = (node.condition_level + indent) * self.config.indent_by
        if indent_by:
            node.whitespaces.value += indent_by

    def visit_ArrayNode(self, node: mparser.ArrayNode) -> None:
        self.enter_node(node)
        if node.args.is_multiline:
            self.level += 1
            if node.args.arguments:
                self.add_nl_after(node.lbracket, indent=self.level)
        node.lbracket.accept(self)
        self.is_function_arguments = False
        node.args.accept(self)
        if node.args.is_multiline:
            self.level -= 1
        node.rbracket.accept(self)
        self.exit_node(node)

    def visit_DictNode(self, node: mparser.DictNode) -> None:
        self.enter_node(node)
        if node.args.is_multiline:
            self.level += 1
            if node.args.kwargs:
                self.add_nl_after(node.lcurl, indent=self.level)
        node.lcurl.accept(self)
        self.is_function_arguments = False
        node.args.accept(self)
        if node.args.is_multiline:
            self.level -= 1
        node.rcurl.accept(self)
        self.exit_node(node)

    def visit_MethodNode(self, node: mparser.MethodNode) -> None:
        self.enter_node(node)
        node.source_object.accept(self)
        is_cont = node.source_object.whitespaces and node.source_object.whitespaces.is_continuation
        if is_cont:
            self.level += 1
        if node.args.is_multiline:
            self.level += 1
            self.add_nl_after(node.lpar, indent=self.level)
        self.is_function_arguments = True
        node.args.accept(self)
        if node.args.is_multiline:
            self.level -= 1
        if is_cont:
            self.level -= 1
        self.exit_node(node)

    def visit_FunctionNode(self, node: mparser.FunctionNode) -> None:
        self.enter_node(node)
        if node.args.is_multiline:
            self.level += 1
            self.add_nl_after(node.lpar, indent=self.level)
        self.is_function_arguments = True
        node.args.accept(self)
        if node.args.is_multiline:
            self.level -= 1
        self.exit_node(node)

    def visit_WhitespaceNode(self, node: mparser.WhitespaceNode) -> None:
        lines = node.value.splitlines(keepends=True)
        if lines:
            indent = (node.condition_level + self.level) * self.config.indent_by
            node.value = '' if node.block_indent else lines.pop(0)
            for line in lines:
                if '#' in line and not line.startswith(indent):
                    node.value += indent
                node.value += line
            if self.indent_after and node.value.endswith(('\n', self.config.indent_by)):
                node.value += indent

    def visit_ArgumentNode(self, node: mparser.ArgumentNode) -> None:
        is_function_arguments = self.is_function_arguments  # record it, because it may change when visiting children
        super().visit_ArgumentNode(node)

        for colon in node.colons:
            self.add_space_after(colon)

        if self.config.wide_colon:
            for key in node.kwargs:
                self.add_space_after(key)

        arguments_count = len(node.arguments) + len(node.kwargs)
        has_trailing_comma = node.commas and len(node.commas) == arguments_count
        if node.is_multiline:
            need_comma = True
            if arguments_count == 1 and is_function_arguments:
                need_comma = not self.config.no_single_comma_function

            if need_comma and not has_trailing_comma:
                comma = mparser.SymbolNode(mparser.Token('comma', node.filename, 0, 0, 0, (0, 0), ','))
                comma.condition_level = node.condition_level
                comma.whitespaces = mparser.WhitespaceNode(mparser.Token('whitespace', node.filename, 0, 0, 0, (0, 0), ''))
                node.commas.append(comma)
            elif has_trailing_comma and not need_comma:
                node.commas.pop(-1)

            arg_index = 0
            if self.config.group_arg_value:
                for arg in node.arguments[:-1]:
                    group_args = False
                    if isinstance(arg, mparser.StringNode) and arg.value.startswith('--'):
                        next_arg = node.arguments[arg_index + 1]
                        if isinstance(next_arg, mparser.StringNode) and not next_arg.value.startswith('--'):
                            group_args = True
                    if group_args:
                        # keep '--arg', 'value' on same line
                        self.add_space_after(node.commas[arg_index])
                    elif arg_index < len(node.commas):
                        self.add_nl_after(node.commas[arg_index], self.level)
                    arg_index += 1

            for comma in node.commas[arg_index:-1]:
                self.add_nl_after(comma, self.level)
            if node.arguments or node.kwargs:
                self.add_nl_after(node, self.level - 1)

        else:
            if has_trailing_comma and not (node.commas[-1].whitespaces and node.commas[-1].whitespaces.value):
                node.commas.pop(-1)

            for comma in node.commas:
                self.add_space_after(comma)

        self.exit_node(node)

    def visit_ParenthesizedNode(self, node: mparser.ParenthesizedNode) -> None:
        self.enter_node(node)
        is_multiline = '\n' in node.lpar.whitespaces.value
        if is_multiline:
            current_indent_after = self.indent_after
            self.indent_after = True
        node.lpar.accept(self)
        node.inner.accept(self)
        if is_multiline:
            self.indent_after = current_indent_after
        node.rpar.accept(self)
        self.exit_node(node)


class ComputeLineLengths(FullAstVisitor):

    def __init__(self, config: FormatterConfig, level: int):
        self.config = config
        self.lengths: T.List[int] = []
        self.length = 0
        self.argument_stack: T.List[mparser.ArgumentNode] = []
        self.level = level
        self.need_regenerate = False

    def visit_default_func(self, node: mparser.BaseNode) -> None:
        self.enter_node(node)
        assert hasattr(node, 'value')
        self.length += len(str(node.value))
        self.exit_node(node)

    def len(self, line: str) -> int:
        '''Compute line length, including tab stops'''
        parts = line.split('\t')
        line_length = len(parts[0])
        for p in parts[1:]:
            tab_length = ((self.length + line_length) % self.config.tab_width) or self.config.tab_width
            line_length += tab_length + len(p)
        return line_length

    def count_multiline(self, value: str) -> None:
        lines = value.splitlines(keepends=True)
        for line in lines:
            if line.endswith('\n'):
                self.lengths.append(self.length + self.len(line) - 1)
                self.length = 0
            else:
                self.length += self.len(line)

    def visit_WhitespaceNode(self, node: mparser.WhitespaceNode) -> None:
        self.count_multiline(node.value)

    def visit_EmptyNode(self, node: mparser.EmptyNode) -> None:
        self.enter_node(node)
        self.exit_node(node)

    def visit_NumberNode(self, node: mparser.NumberNode) -> None:
        self.enter_node(node)
        self.length += len(node.raw_value)
        self.exit_node(node)

    def visit_StringNode(self, node: mparser.StringNode) -> None:
        self.enter_node(node)
        if node.is_fstring:
            self.length += 1

        if node.is_multiline:
            self.length += 3
            self.count_multiline(node.value)
            self.length += 3
        else:
            self.length += self.len(node.raw_value) + 2

        self.exit_node(node)

    def visit_ContinueNode(self, node: mparser.ContinueNode) -> None:
        self.enter_node(node)
        self.length += len('continue')
        self.exit_node(node)

    def visit_BreakNode(self, node: mparser.BreakNode) -> None:
        self.enter_node(node)
        self.length += len('break')
        self.exit_node(node)

    def split_if_needed(self, node: mparser.ArgumentNode) -> None:
        if len(node) and not node.is_multiline and self.length > self.config.max_line_length:
            arg = self.argument_stack[self.level] if len(self.argument_stack) > self.level else node
            arg.is_multiline = True
            self.need_regenerate = True

    def visit_ArgumentNode(self, node: mparser.ArgumentNode) -> None:
        self.argument_stack.append(node)
        super().visit_ArgumentNode(node)
        self.split_if_needed(node)
        self.argument_stack.pop(-1)

    def visit_ArrayNode(self, node: mparser.ArrayNode) -> None:
        self.enter_node(node)
        node.lbracket.accept(self)
        node.args.accept(self)
        node.rbracket.accept(self)
        self.split_if_needed(node.args)  # split if closing bracket is too far
        self.exit_node(node)

    def visit_DictNode(self, node: mparser.DictNode) -> None:
        self.enter_node(node)
        node.lcurl.accept(self)
        node.args.accept(self)
        node.rcurl.accept(self)
        self.split_if_needed(node.args)  # split if closing bracket is too far
        self.exit_node(node)


class SubdirFetcher(FullAstVisitor):

    def __init__(self, current_dir: Path):
        self.current_dir = current_dir
        self.subdirs: T.List[Path] = []

    def visit_FunctionNode(self, node: mparser.FunctionNode) -> None:
        if node.func_name.value == 'subdir':
            if node.args.arguments and isinstance(node.args.arguments[0], mparser.StringNode):
                subdir = node.args.arguments[0].value
                self.subdirs.append(self.current_dir / subdir)
        super().visit_FunctionNode(node)


class Formatter:

    def __init__(self, configuration_file: T.Optional[Path], use_editor_config: bool, fetch_subdirs: bool):
        self.fetch_subdirs = fetch_subdirs
        self.use_editor_config = use_editor_config
        self.config = self.load_configuration(configuration_file)
        self.current_config = self.config

        self.current_dir = Path()
        self.subdirs: T.List[Path] = []

    def load_editor_config(self, source_file: Path) -> EditorConfig:
        # See https://editorconfig.org/
        config = EditorConfig()

        for p in source_file.parents:
            editorconfig_file = p / '.editorconfig'
            if not editorconfig_file.exists():
                continue

            cp = DefaultConfigParser(delimiters=('=',))
            cp.read_default(editorconfig_file)

            sections = [section for section in cp.sections() if match_path(source_file.as_posix(), section)]
            for f in fields(config):
                if getattr(cp, f.name, None) is not None:
                    continue  # value already set from higher file

                getter = f.metadata['getter']
                for section in sections:
                    try:
                        value = getter(cp, section, f.name, fallback=None)
                    except ValueError as e:
                        raise MesonException(f'Invalid type for key "{f.name}" in "{editorconfig_file}" file:\n{e}') from e
                    if value is not None:
                        setattr(config, f.name, value)

            # Root is not required except in the top level .editorconfig.
            if cp.getboolean(cp.default_section, 'root', fallback=False):
                break

        return config

    def load_configuration(self, configuration_file: T.Optional[Path]) -> FormatterConfig:
        config = FormatterConfig()
        if configuration_file:
            cp = DefaultConfigParser()
            try:
                cp.read_default(configuration_file)
            except ParsingError as e:
                raise MesonException(f'Unable to parse configuration file "{configuration_file}":\n{e}') from e

            extra_keys = sorted(set(cp.defaults()).difference(f.name for f in fields(config)))
            if extra_keys:
                raise MesonException(f'Unknown config keys: "{", ".join(extra_keys)}" in configuration file "{configuration_file}"')

            for f in fields(config):
                getter = f.metadata['getter']
                try:
                    value = getter(cp, cp.default_section, f.name, fallback=None)
                except ValueError as e:
                    raise MesonException(
                        f'Error parsing "{str(configuration_file)}", option "{f.name}", error: "{e!s}"')
                if value is not None:
                    setattr(config, f.name, value)

            if config.use_editor_config:
                self.use_editor_config = True

        return config

    def format(self, code: str, source_file: Path) -> str:
        self.current_dir = source_file.parent
        self.current_config = FormatterConfig.default()
        if self.use_editor_config:
            self.current_config = self.current_config.with_editorconfig(self.load_editor_config(source_file))
        self.current_config = self.current_config.update(self.config)

        ast = mparser.Parser(code, source_file.as_posix()).parse()
        if self.fetch_subdirs:
            subdir_fetcher = SubdirFetcher(self.current_dir)
            ast.accept(subdir_fetcher)
            self.subdirs = subdir_fetcher.subdirs

        ast.accept(AstConditionLevel())
        for level in range(5):
            ast.accept(TrimWhitespaces(self.current_config))
            ast.accept(ArgumentFormatter(self.current_config))

            cll = ComputeLineLengths(self.current_config, level)
            ast.accept(cll)
            if not cll.need_regenerate:
                break

        printer = RawPrinter()
        ast.accept(printer)
        return printer.result


def add_arguments(parser: argparse.ArgumentParser) -> None:
    inplace_group = parser.add_mutually_exclusive_group()
    inplace_group.add_argument(
        '-q', '--check-only',
        action='store_true',
        help='exit with 1 if files would be modified by meson format'
    )
    inplace_group.add_argument(
        '-i', '--inplace',
        action='store_true',
        help='format files in-place'
    )
    parser.add_argument(
        '-r', '--recursive',
        action='store_true',
        help='recurse subdirs (requires --check-only or --inplace option)',
    )
    parser.add_argument(
        '-c', '--configuration',
        metavar='meson.format',
        type=Path,
        help='read configuration from meson.format'
    )
    parser.add_argument(
        '-e', '--editor-config',
        action='store_true',
        default=False,
        help='try to read configuration from .editorconfig'
    )
    parser.add_argument(
        '-o', '--output',
        type=Path,
        help='output file (implies having exactly one input)'
    )
    parser.add_argument(
        'sources',
        nargs='*',
        type=Path,
        help='meson source files'
    )

def run(options: argparse.Namespace) -> int:
    if options.output and len(options.sources) != 1:
        raise MesonException('--output argument implies having exactly one source file')
    if options.recursive and not (options.inplace or options.check_only):
        raise MesonException('--recursive argument requires either --inplace or --check-only option')

    from_stdin = len(options.sources) == 1 and options.sources[0].name == '-' and options.sources[0].parent == Path()
    if options.recursive and from_stdin:
        raise MesonException('--recursive argument is not compatible with stdin input')
    if options.inplace and from_stdin:
        raise MesonException('--inplace argument is not compatible with stdin input')

    sources: T.List[Path] = options.sources.copy() or [Path(build_filename)]
    if not options.configuration:
        default_config_path = sources[0].parent / 'meson.format'
        if default_config_path.exists():
            options.configuration = default_config_path
    formatter = Formatter(options.configuration, options.editor_config, options.recursive)

    while sources:
        src_file = sources.pop(0)
        if src_file.is_dir():
            src_file = src_file / build_filename

        try:
            if from_stdin:
                src_file = Path('STDIN')  # used for error messages and introspection
                code = sys.stdin.read()
            else:
                code = src_file.read_text(encoding='utf-8')
        except IOError as e:
            raise MesonException(f'Unable to read from {src_file}') from e

        formatted = formatter.format(code, src_file)
        if options.recursive:
            sources.extend(formatter.subdirs)

        if options.inplace:
            try:
                with src_file.open('w', encoding='utf-8', newline=formatter.current_config.newline) as sf:
                    sf.write(formatted)
            except IOError as e:
                raise MesonException(f'Unable to write to {src_file}') from e
        elif options.check_only:
            # TODO: add verbose output showing diffs
            if code != formatted:
                return 1
        elif options.output:
            try:
                with options.output.open('w', encoding='utf-8', newline=formatter.current_config.newline) as of:
                    of.write(formatted)
            except IOError as e:
                raise MesonException(f'Unable to write to {options.output}') from e
        else:
            print(formatted, end='')

    return 0

# TODO: remove empty newlines when more than N (2...)
# TODO: magic comment to prevent formatting
# TODO: split long lines on binary operators
# TODO: align series of assignments
# TODO: align comments
# TODO: move comments on long lines

# Differences from muon format:
# - By default, uses two spaces before comment, and added option for that
# - Muon will mix CRLF and LF on Windows files...
# - Support for end_of_line char
# - Support for max_line_length, end_of_line, insert_final_newline, tab_width in .editorconfig
# - Option to simplify string literals
# - Option to recognize and parse meson.build in subdirs
# - Correctly compute line length when using tabs
# - By default, arguments in files() are sorted alphabetically
# - Option to group '--arg', 'value' on same line in multiline arguments
