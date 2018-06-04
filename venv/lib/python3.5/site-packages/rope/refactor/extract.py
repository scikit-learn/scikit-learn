import re

from rope.base.utils.datastructures import OrderedSet
from rope.base import ast, codeanalyze
from rope.base.change import ChangeSet, ChangeContents
from rope.base.exceptions import RefactoringError
from rope.base.utils import pycompat
from rope.refactor import (sourceutils, similarfinder,
                           patchedast, suites, usefunction)


# Extract refactoring has lots of special cases.  I tried to split it
# to smaller parts to make it more manageable:
#
# _ExtractInfo: holds information about the refactoring; it is passed
# to the parts that need to have information about the refactoring
#
# _ExtractCollector: merely saves all of the information necessary for
# performing the refactoring.
#
# _DefinitionLocationFinder: finds where to insert the definition.
#
# _ExceptionalConditionChecker: checks for exceptional conditions in
# which the refactoring cannot be applied.
#
# _ExtractMethodParts: generates the pieces of code (like definition)
# needed for performing extract method.
#
# _ExtractVariableParts: like _ExtractMethodParts for variables.
#
# _ExtractPerformer: Uses above classes to collect refactoring
# changes.
#
# There are a few more helper functions and classes used by above
# classes.
class _ExtractRefactoring(object):

    def __init__(self, project, resource, start_offset, end_offset,
                 variable=False):
        self.project = project
        self.resource = resource
        self.start_offset = self._fix_start(resource.read(), start_offset)
        self.end_offset = self._fix_end(resource.read(), end_offset)

    def _fix_start(self, source, offset):
        while offset < len(source) and source[offset].isspace():
            offset += 1
        return offset

    def _fix_end(self, source, offset):
        while offset > 0 and source[offset - 1].isspace():
            offset -= 1
        return offset

    def get_changes(self, extracted_name, similar=False, global_=False):
        """Get the changes this refactoring makes

        :parameters:
            - `similar`: if `True`, similar expressions/statements are also
              replaced.
            - `global_`: if `True`, the extracted method/variable will
              be global.

        """
        info = _ExtractInfo(
            self.project, self.resource, self.start_offset, self.end_offset,
            extracted_name, variable=self.kind == 'variable',
            similar=similar, make_global=global_)
        new_contents = _ExtractPerformer(info).extract()
        changes = ChangeSet('Extract %s <%s>' % (self.kind,
                                                 extracted_name))
        changes.add_change(ChangeContents(self.resource, new_contents))
        return changes


class ExtractMethod(_ExtractRefactoring):

    def __init__(self, *args, **kwds):
        super(ExtractMethod, self).__init__(*args, **kwds)

    kind = 'method'


class ExtractVariable(_ExtractRefactoring):

    def __init__(self, *args, **kwds):
        kwds = dict(kwds)
        kwds['variable'] = True
        super(ExtractVariable, self).__init__(*args, **kwds)

    kind = 'variable'


class _ExtractInfo(object):
    """Holds information about the extract to be performed"""

    def __init__(self, project, resource, start, end, new_name,
                 variable, similar, make_global):
        self.project = project
        self.resource = resource
        self.pymodule = project.get_pymodule(resource)
        self.global_scope = self.pymodule.get_scope()
        self.source = self.pymodule.source_code
        self.lines = self.pymodule.lines
        self.new_name = new_name
        self.variable = variable
        self.similar = similar
        self._init_parts(start, end)
        self._init_scope()
        self.make_global = make_global

    def _init_parts(self, start, end):
        self.region = (self._choose_closest_line_end(start),
                       self._choose_closest_line_end(end, end=True))

        start = self.logical_lines.logical_line_in(
            self.lines.get_line_number(self.region[0]))[0]
        end = self.logical_lines.logical_line_in(
            self.lines.get_line_number(self.region[1]))[1]
        self.region_lines = (start, end)

        self.lines_region = (self.lines.get_line_start(self.region_lines[0]),
                             self.lines.get_line_end(self.region_lines[1]))

    @property
    def logical_lines(self):
        return self.pymodule.logical_lines

    def _init_scope(self):
        start_line = self.region_lines[0]
        scope = self.global_scope.get_inner_scope_for_line(start_line)
        if scope.get_kind() != 'Module' and scope.get_start() == start_line:
            scope = scope.parent
        self.scope = scope
        self.scope_region = self._get_scope_region(self.scope)

    def _get_scope_region(self, scope):
        return (self.lines.get_line_start(scope.get_start()),
                self.lines.get_line_end(scope.get_end()) + 1)

    def _choose_closest_line_end(self, offset, end=False):
        lineno = self.lines.get_line_number(offset)
        line_start = self.lines.get_line_start(lineno)
        line_end = self.lines.get_line_end(lineno)
        if self.source[line_start:offset].strip() == '':
            if end:
                return line_start - 1
            else:
                return line_start
        elif self.source[offset:line_end].strip() == '':
            return min(line_end, len(self.source))
        return offset

    @property
    def one_line(self):
        return self.region != self.lines_region and \
            (self.logical_lines.logical_line_in(self.region_lines[0]) ==
             self.logical_lines.logical_line_in(self.region_lines[1]))

    @property
    def global_(self):
        return self.scope.parent is None

    @property
    def method(self):
        return self.scope.parent is not None and \
            self.scope.parent.get_kind() == 'Class'

    @property
    def indents(self):
        return sourceutils.get_indents(self.pymodule.lines,
                                       self.region_lines[0])

    @property
    def scope_indents(self):
        if self.global_:
            return 0
        return sourceutils.get_indents(self.pymodule.lines,
                                       self.scope.get_start())

    @property
    def extracted(self):
        return self.source[self.region[0]:self.region[1]]

    _returned = None

    @property
    def returned(self):
        """Does the extracted piece contain return statement"""
        if self._returned is None:
            node = _parse_text(self.extracted)
            self._returned = usefunction._returns_last(node)
        return self._returned


class _ExtractCollector(object):
    """Collects information needed for performing the extract"""

    def __init__(self, info):
        self.definition = None
        self.body_pattern = None
        self.checks = {}
        self.replacement_pattern = None
        self.matches = None
        self.replacements = None
        self.definition_location = None


class _ExtractPerformer(object):

    def __init__(self, info):
        self.info = info
        _ExceptionalConditionChecker()(self.info)

    def extract(self):
        extract_info = self._collect_info()
        content = codeanalyze.ChangeCollector(self.info.source)
        definition = extract_info.definition
        lineno, indents = extract_info.definition_location
        offset = self.info.lines.get_line_start(lineno)
        indented = sourceutils.fix_indentation(definition, indents)
        content.add_change(offset, offset, indented)
        self._replace_occurrences(content, extract_info)
        return content.get_changed()

    def _replace_occurrences(self, content, extract_info):
        for match in extract_info.matches:
            replacement = similarfinder.CodeTemplate(
                extract_info.replacement_pattern)
            mapping = {}
            for name in replacement.get_names():
                node = match.get_ast(name)
                if node:
                    start, end = patchedast.node_region(match.get_ast(name))
                    mapping[name] = self.info.source[start:end]
                else:
                    mapping[name] = name
            region = match.get_region()
            content.add_change(region[0], region[1],
                               replacement.substitute(mapping))

    def _collect_info(self):
        extract_collector = _ExtractCollector(self.info)
        self._find_definition(extract_collector)
        self._find_matches(extract_collector)
        self._find_definition_location(extract_collector)
        return extract_collector

    def _find_matches(self, collector):
        regions = self._where_to_search()
        finder = similarfinder.SimilarFinder(self.info.pymodule)
        matches = []
        for start, end in regions:
            matches.extend((finder.get_matches(collector.body_pattern,
                                               collector.checks, start, end)))
        collector.matches = matches

    def _where_to_search(self):
        if self.info.similar:
            if self.info.make_global or self.info.global_:
                return [(0, len(self.info.pymodule.source_code))]
            if self.info.method and not self.info.variable:
                class_scope = self.info.scope.parent
                regions = []
                method_kind = _get_function_kind(self.info.scope)
                for scope in class_scope.get_scopes():
                    if method_kind == 'method' and \
                       _get_function_kind(scope) != 'method':
                        continue
                    start = self.info.lines.get_line_start(scope.get_start())
                    end = self.info.lines.get_line_end(scope.get_end())
                    regions.append((start, end))
                return regions
            else:
                if self.info.variable:
                    return [self.info.scope_region]
                else:
                    return [self.info._get_scope_region(
                        self.info.scope.parent)]
        else:
            return [self.info.region]

    def _find_definition_location(self, collector):
        matched_lines = []
        for match in collector.matches:
            start = self.info.lines.get_line_number(match.get_region()[0])
            start_line = self.info.logical_lines.logical_line_in(start)[0]
            matched_lines.append(start_line)
        location_finder = _DefinitionLocationFinder(self.info, matched_lines)
        collector.definition_location = (location_finder.find_lineno(),
                                         location_finder.find_indents())

    def _find_definition(self, collector):
        if self.info.variable:
            parts = _ExtractVariableParts(self.info)
        else:
            parts = _ExtractMethodParts(self.info)
        collector.definition = parts.get_definition()
        collector.body_pattern = parts.get_body_pattern()
        collector.replacement_pattern = parts.get_replacement_pattern()
        collector.checks = parts.get_checks()


class _DefinitionLocationFinder(object):

    def __init__(self, info, matched_lines):
        self.info = info
        self.matched_lines = matched_lines
        # This only happens when subexpressions cannot be matched
        if not matched_lines:
            self.matched_lines.append(self.info.region_lines[0])

    def find_lineno(self):
        if self.info.variable and not self.info.make_global:
            return self._get_before_line()
        if self.info.make_global or self.info.global_:
            toplevel = self._find_toplevel(self.info.scope)
            ast = self.info.pymodule.get_ast()
            newlines = sorted(self.matched_lines + [toplevel.get_end() + 1])
            return suites.find_visible(ast, newlines)
        return self._get_after_scope()

    def _find_toplevel(self, scope):
        toplevel = scope
        if toplevel.parent is not None:
            while toplevel.parent.parent is not None:
                toplevel = toplevel.parent
        return toplevel

    def find_indents(self):
        if self.info.variable and not self.info.make_global:
            return sourceutils.get_indents(self.info.lines,
                                           self._get_before_line())
        else:
            if self.info.global_ or self.info.make_global:
                return 0
        return self.info.scope_indents

    def _get_before_line(self):
        ast = self.info.scope.pyobject.get_ast()
        return suites.find_visible(ast, self.matched_lines)

    def _get_after_scope(self):
        return self.info.scope.get_end() + 1


class _ExceptionalConditionChecker(object):

    def __call__(self, info):
        self.base_conditions(info)
        if info.one_line:
            self.one_line_conditions(info)
        else:
            self.multi_line_conditions(info)

    def base_conditions(self, info):
        if info.region[1] > info.scope_region[1]:
            raise RefactoringError('Bad region selected for extract method')
        end_line = info.region_lines[1]
        end_scope = info.global_scope.get_inner_scope_for_line(end_line)
        if end_scope != info.scope and end_scope.get_end() != end_line:
            raise RefactoringError('Bad region selected for extract method')
        try:
            extracted = info.source[info.region[0]:info.region[1]]
            if info.one_line:
                extracted = '(%s)' % extracted
            if _UnmatchedBreakOrContinueFinder.has_errors(extracted):
                raise RefactoringError('A break/continue without having a '
                                       'matching for/while loop.')
        except SyntaxError:
            raise RefactoringError('Extracted piece should '
                                   'contain complete statements.')

    def one_line_conditions(self, info):
        if self._is_region_on_a_word(info):
            raise RefactoringError('Should extract complete statements.')
        if info.variable and not info.one_line:
            raise RefactoringError('Extract variable should not '
                                   'span multiple lines.')

    def multi_line_conditions(self, info):
        node = _parse_text(info.source[info.region[0]:info.region[1]])
        count = usefunction._return_count(node)
        if count > 1:
            raise RefactoringError('Extracted piece can have only one '
                                   'return statement.')
        if usefunction._yield_count(node):
            raise RefactoringError('Extracted piece cannot '
                                   'have yield statements.')
        if count == 1 and not usefunction._returns_last(node):
            raise RefactoringError('Return should be the last statement.')
        if info.region != info.lines_region:
            raise RefactoringError('Extracted piece should '
                                   'contain complete statements.')

    def _is_region_on_a_word(self, info):
        if info.region[0] > 0 and \
                self._is_on_a_word(info, info.region[0] - 1) or \
                self._is_on_a_word(info, info.region[1] - 1):
            return True

    def _is_on_a_word(self, info, offset):
        prev = info.source[offset]
        if not (prev.isalnum() or prev == '_') or \
           offset + 1 == len(info.source):
            return False
        next = info.source[offset + 1]
        return next.isalnum() or next == '_'


class _ExtractMethodParts(object):

    def __init__(self, info):
        self.info = info
        self.info_collector = self._create_info_collector()

    def get_definition(self):
        if self.info.global_:
            return '\n%s\n' % self._get_function_definition()
        else:
            return '\n%s' % self._get_function_definition()

    def get_replacement_pattern(self):
        variables = []
        variables.extend(self._find_function_arguments())
        variables.extend(self._find_function_returns())
        return similarfinder.make_pattern(self._get_call(), variables)

    def get_body_pattern(self):
        variables = []
        variables.extend(self._find_function_arguments())
        variables.extend(self._find_function_returns())
        variables.extend(self._find_temps())
        return similarfinder.make_pattern(self._get_body(), variables)

    def _get_body(self):
        result = sourceutils.fix_indentation(self.info.extracted, 0)
        if self.info.one_line:
            result = '(%s)' % result
        return result

    def _find_temps(self):
        return usefunction.find_temps(self.info.project,
                                      self._get_body())

    def get_checks(self):
        if self.info.method and not self.info.make_global:
            if _get_function_kind(self.info.scope) == 'method':
                class_name = similarfinder._pydefined_to_str(
                    self.info.scope.parent.pyobject)
                return {self._get_self_name(): 'type=' + class_name}
        return {}

    def _create_info_collector(self):
        zero = self.info.scope.get_start() - 1
        start_line = self.info.region_lines[0] - zero
        end_line = self.info.region_lines[1] - zero
        info_collector = _FunctionInformationCollector(start_line, end_line,
                                                       self.info.global_)
        body = self.info.source[self.info.scope_region[0]:
                                self.info.scope_region[1]]
        node = _parse_text(body)
        ast.walk(node, info_collector)
        return info_collector

    def _get_function_definition(self):
        args = self._find_function_arguments()
        returns = self._find_function_returns()
        result = []
        if self.info.method and not self.info.make_global and \
           _get_function_kind(self.info.scope) != 'method':
            result.append('@staticmethod\n')
        result.append('def %s:\n' % self._get_function_signature(args))
        unindented_body = self._get_unindented_function_body(returns)
        indents = sourceutils.get_indent(self.info.project)
        function_body = sourceutils.indent_lines(unindented_body, indents)
        result.append(function_body)
        definition = ''.join(result)

        return definition + '\n'

    def _get_function_signature(self, args):
        args = list(args)
        prefix = ''
        if self._extracting_method():
            self_name = self._get_self_name()
            if self_name is None:
                raise RefactoringError('Extracting a method from a function '
                                       'with no self argument.')
            if self_name in args:
                args.remove(self_name)
            args.insert(0, self_name)
        return prefix + self.info.new_name + \
            '(%s)' % self._get_comma_form(args)

    def _extracting_method(self):
        return self.info.method and not self.info.make_global and \
            _get_function_kind(self.info.scope) == 'method'

    def _get_self_name(self):
        param_names = self.info.scope.pyobject.get_param_names()
        if param_names:
            return param_names[0]

    def _get_function_call(self, args):
        prefix = ''
        if self.info.method and not self.info.make_global:
            if _get_function_kind(self.info.scope) == 'method':
                self_name = self._get_self_name()
                if self_name in args:
                    args.remove(self_name)
                prefix = self_name + '.'
            else:
                prefix = self.info.scope.parent.pyobject.get_name() + '.'
        return prefix + '%s(%s)' % (self.info.new_name,
                                    self._get_comma_form(args))

    def _get_comma_form(self, names):
        result = ''
        if names:
            result += names[0]
            for name in names[1:]:
                result += ', ' + name
        return result

    def _get_call(self):
        if self.info.one_line:
            args = self._find_function_arguments()
            return self._get_function_call(args)
        args = self._find_function_arguments()
        returns = self._find_function_returns()
        call_prefix = ''
        if returns:
            call_prefix = self._get_comma_form(returns) + ' = '
        if self.info.returned:
            call_prefix = 'return '
        return call_prefix + self._get_function_call(args)

    def _find_function_arguments(self):
        # if not make_global, do not pass any global names; they are
        # all visible.
        if self.info.global_ and not self.info.make_global:
            return ()
        if not self.info.one_line:
            result = (self.info_collector.prewritten &
                      self.info_collector.read)
            result |= (self.info_collector.prewritten &
                       self.info_collector.postread &
                       (self.info_collector.maybe_written -
                        self.info_collector.written))
            return list(result)
        start = self.info.region[0]
        if start == self.info.lines_region[0]:
            start = start + re.search('\S', self.info.extracted).start()
        function_definition = self.info.source[start:self.info.region[1]]
        read = _VariableReadsAndWritesFinder.find_reads_for_one_liners(
            function_definition)
        return list(self.info_collector.prewritten.intersection(read))

    def _find_function_returns(self):
        if self.info.one_line or self.info.returned:
            return []
        written = self.info_collector.written | \
            self.info_collector.maybe_written
        return list(written & self.info_collector.postread)

    def _get_unindented_function_body(self, returns):
        if self.info.one_line:
            return 'return ' + _join_lines(self.info.extracted)
        extracted_body = self.info.extracted
        unindented_body = sourceutils.fix_indentation(extracted_body, 0)
        if returns:
            unindented_body += '\nreturn %s' % self._get_comma_form(returns)
        return unindented_body


class _ExtractVariableParts(object):

    def __init__(self, info):
        self.info = info

    def get_definition(self):
        result = self.info.new_name + ' = ' + \
            _join_lines(self.info.extracted) + '\n'
        return result

    def get_body_pattern(self):
        return '(%s)' % self.info.extracted.strip()

    def get_replacement_pattern(self):
        return self.info.new_name

    def get_checks(self):
        return {}


class _FunctionInformationCollector(object):

    def __init__(self, start, end, is_global):
        self.start = start
        self.end = end
        self.is_global = is_global
        self.prewritten = OrderedSet()
        self.maybe_written = OrderedSet()
        self.written = OrderedSet()
        self.read = OrderedSet()
        self.postread = OrderedSet()
        self.postwritten = OrderedSet()
        self.host_function = True
        self.conditional = False

    def _read_variable(self, name, lineno):
        if self.start <= lineno <= self.end:
            if name not in self.written:
                if not self.conditional or name not in self.maybe_written:
                    self.read.add(name)
        if self.end < lineno:
            if name not in self.postwritten:
                self.postread.add(name)

    def _written_variable(self, name, lineno):
        if self.start <= lineno <= self.end:
            if self.conditional:
                self.maybe_written.add(name)
            else:
                self.written.add(name)
        if self.start > lineno:
            self.prewritten.add(name)
        if self.end < lineno:
            self.postwritten.add(name)

    def _FunctionDef(self, node):
        if not self.is_global and self.host_function:
            self.host_function = False
            for name in _get_argnames(node.args):
                self._written_variable(name, node.lineno)
            for child in node.body:
                ast.walk(child, self)
        else:
            self._written_variable(node.name, node.lineno)
            visitor = _VariableReadsAndWritesFinder()
            for child in node.body:
                ast.walk(child, visitor)
            for name in visitor.read - visitor.written:
                self._read_variable(name, node.lineno)

    def _Name(self, node):
        if isinstance(node.ctx, (ast.Store, ast.AugStore)):
            self._written_variable(node.id, node.lineno)
        if not isinstance(node.ctx, ast.Store):
            self._read_variable(node.id, node.lineno)

    def _Assign(self, node):
        ast.walk(node.value, self)
        for child in node.targets:
            ast.walk(child, self)

    def _ClassDef(self, node):
        self._written_variable(node.name, node.lineno)

    def _handle_conditional_node(self, node):
        self.conditional = True
        try:
            for child in ast.get_child_nodes(node):
                ast.walk(child, self)
        finally:
            self.conditional = False

    def _If(self, node):
        self._handle_conditional_node(node)

    def _While(self, node):
        self._handle_conditional_node(node)

    def _For(self, node):
        self.conditional = True
        try:
            # iter has to be checked before the target variables
            ast.walk(node.iter, self)
            ast.walk(node.target, self)

            for child in node.body:
                ast.walk(child, self)
            for child in node.orelse:
                ast.walk(child, self)
        finally:
            self.conditional = False


def _get_argnames(arguments):
    result = [pycompat.get_ast_arg_arg(node) for node in arguments.args
              if isinstance(node, pycompat.ast_arg_type)]
    if arguments.vararg:
        result.append(pycompat.get_ast_arg_arg(arguments.vararg))
    if arguments.kwarg:
        result.append(pycompat.get_ast_arg_arg(arguments.kwarg))
    return result


class _VariableReadsAndWritesFinder(object):

    def __init__(self):
        self.written = set()
        self.read = set()

    def _Name(self, node):
        if isinstance(node.ctx, (ast.Store, ast.AugStore)):
            self.written.add(node.id)
        if not isinstance(node, ast.Store):
            self.read.add(node.id)

    def _FunctionDef(self, node):
        self.written.add(node.name)
        visitor = _VariableReadsAndWritesFinder()
        for child in ast.get_child_nodes(node):
            ast.walk(child, visitor)
        self.read.update(visitor.read - visitor.written)

    def _Class(self, node):
        self.written.add(node.name)

    @staticmethod
    def find_reads_and_writes(code):
        if code.strip() == '':
            return set(), set()
        if isinstance(code, unicode):
            code = code.encode('utf-8')
        node = _parse_text(code)
        visitor = _VariableReadsAndWritesFinder()
        ast.walk(node, visitor)
        return visitor.read, visitor.written

    @staticmethod
    def find_reads_for_one_liners(code):
        if code.strip() == '':
            return set(), set()
        node = _parse_text(code)
        visitor = _VariableReadsAndWritesFinder()
        ast.walk(node, visitor)
        return visitor.read


class _UnmatchedBreakOrContinueFinder(object):

    def __init__(self):
        self.error = False
        self.loop_count = 0

    def _For(self, node):
        self.loop_encountered(node)

    def _While(self, node):
        self.loop_encountered(node)

    def loop_encountered(self, node):
        self.loop_count += 1
        for child in node.body:
            ast.walk(child, self)
        self.loop_count -= 1
        if node.orelse:
            if isinstance(node.orelse,(list,tuple)):
                for node_ in node.orelse:
                    ast.walk(node_, self)
            else:
                ast.walk(node.orelse, self)

    def _Break(self, node):
        self.check_loop()

    def _Continue(self, node):
        self.check_loop()

    def check_loop(self):
        if self.loop_count < 1:
            self.error = True

    def _FunctionDef(self, node):
        pass

    def _ClassDef(self, node):
        pass

    @staticmethod
    def has_errors(code):
        if code.strip() == '':
            return False
        node = _parse_text(code)
        visitor = _UnmatchedBreakOrContinueFinder()
        ast.walk(node, visitor)
        return visitor.error


def _get_function_kind(scope):
    return scope.pyobject.get_kind()


def _parse_text(body):
    body = sourceutils.fix_indentation(body, 0)
    node = ast.parse(body)
    return node


def _join_lines(code):
    lines = []
    for line in code.splitlines():
        if line.endswith('\\'):
            lines.append(line[:-1].strip())
        else:
            lines.append(line.strip())
    return ' '.join(lines)
