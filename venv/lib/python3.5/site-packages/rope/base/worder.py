import bisect
import keyword

import rope.base.simplify


def get_name_at(resource, offset):
    source_code = resource.read()
    word_finder = Worder(source_code)
    return word_finder.get_word_at(offset)


class Worder(object):
    """A class for finding boundaries of words and expressions

    Note that in these methods, offset should be the index of the
    character not the index of the character after it.
    """

    def __init__(self, code, handle_ignores=False):
        simplified = rope.base.simplify.real_code(code)
        self.code_finder = _RealFinder(simplified, code)
        self.handle_ignores = handle_ignores
        self.code = code

    def _init_ignores(self):
        ignores = rope.base.simplify.ignored_regions(self.code)
        self.dumb_finder = _RealFinder(self.code, self.code)
        self.starts = [ignored[0] for ignored in ignores]
        self.ends = [ignored[1] for ignored in ignores]

    def _context_call(self, name, offset):
        if self.handle_ignores:
            if not hasattr(self, 'starts'):
                self._init_ignores()
            start = bisect.bisect(self.starts, offset)
            if start > 0 and offset < self.ends[start - 1]:
                return getattr(self.dumb_finder, name)(offset)
        return getattr(self.code_finder, name)(offset)

    def get_primary_at(self, offset):
        return self._context_call('get_primary_at', offset)

    def get_word_at(self, offset):
        return self._context_call('get_word_at', offset)

    def get_primary_range(self, offset):
        return self._context_call('get_primary_range', offset)

    def get_splitted_primary_before(self, offset):
        return self._context_call('get_splitted_primary_before', offset)

    def get_word_range(self, offset):
        return self._context_call('get_word_range', offset)

    def is_function_keyword_parameter(self, offset):
        return self.code_finder.is_function_keyword_parameter(offset)

    def is_a_class_or_function_name_in_header(self, offset):
        return self.code_finder.is_a_class_or_function_name_in_header(offset)

    def is_from_statement_module(self, offset):
        return self.code_finder.is_from_statement_module(offset)

    def is_from_aliased(self, offset):
        return self.code_finder.is_from_aliased(offset)

    def find_parens_start_from_inside(self, offset):
        return self.code_finder.find_parens_start_from_inside(offset)

    def is_a_name_after_from_import(self, offset):
        return self.code_finder.is_a_name_after_from_import(offset)

    def is_from_statement(self, offset):
        return self.code_finder.is_from_statement(offset)

    def get_from_aliased(self, offset):
        return self.code_finder.get_from_aliased(offset)

    def is_import_statement(self, offset):
        return self.code_finder.is_import_statement(offset)

    def is_assigned_here(self, offset):
        return self.code_finder.is_assigned_here(offset)

    def is_a_function_being_called(self, offset):
        return self.code_finder.is_a_function_being_called(offset)

    def get_word_parens_range(self, offset):
        return self.code_finder.get_word_parens_range(offset)

    def is_name_assigned_in_class_body(self, offset):
        return self.code_finder.is_name_assigned_in_class_body(offset)

    def is_on_function_call_keyword(self, offset):
        return self.code_finder.is_on_function_call_keyword(offset)

    def _find_parens_start(self, offset):
        return self.code_finder._find_parens_start(offset)

    def get_parameters(self, first, last):
        return self.code_finder.get_parameters(first, last)

    def get_from_module(self, offset):
        return self.code_finder.get_from_module(offset)

    def is_assigned_in_a_tuple_assignment(self, offset):
        return self.code_finder.is_assigned_in_a_tuple_assignment(offset)

    def get_assignment_type(self, offset):
        return self.code_finder.get_assignment_type(offset)

    def get_function_and_args_in_header(self, offset):
        return self.code_finder.get_function_and_args_in_header(offset)

    def get_lambda_and_args(self, offset):
        return self.code_finder.get_lambda_and_args(offset)

    def find_function_offset(self, offset):
        return self.code_finder.find_function_offset(offset)


class _RealFinder(object):

    def __init__(self, code, raw):
        self.code = code
        self.raw = raw

    def _find_word_start(self, offset):
        current_offset = offset
        while current_offset >= 0 and self._is_id_char(current_offset):
            current_offset -= 1
        return current_offset + 1

    def _find_word_end(self, offset):
        while offset + 1 < len(self.code) and self._is_id_char(offset + 1):
            offset += 1
        return offset

    def _find_last_non_space_char(self, offset):
        while offset >= 0 and self.code[offset].isspace():
            if self.code[offset] == '\n':
                return offset
            offset -= 1
        return max(-1, offset)

    def get_word_at(self, offset):
        offset = self._get_fixed_offset(offset)
        return self.raw[self._find_word_start(offset):
                        self._find_word_end(offset) + 1]

    def _get_fixed_offset(self, offset):
        if offset >= len(self.code):
            return offset - 1
        if not self._is_id_char(offset):
            if offset > 0 and self._is_id_char(offset - 1):
                return offset - 1
            if offset < len(self.code) - 1 and self._is_id_char(offset + 1):
                return offset + 1
        return offset

    def _is_id_char(self, offset):
        return self.code[offset].isalnum() or self.code[offset] == '_'

    def _find_string_start(self, offset):
        kind = self.code[offset]
        try:
            return self.code.rindex(kind, 0, offset)
        except ValueError:
            return 0

    def _find_parens_start(self, offset):
        offset = self._find_last_non_space_char(offset - 1)
        while offset >= 0 and self.code[offset] not in '[({':
            if self.code[offset] not in ':,':
                offset = self._find_primary_start(offset)
            offset = self._find_last_non_space_char(offset - 1)
        return offset

    def _find_atom_start(self, offset):
        old_offset = offset
        if self.code[offset] == '\n':
            return offset + 1
        if self.code[offset].isspace():
            offset = self._find_last_non_space_char(offset)
        if self.code[offset] in '\'"':
            return self._find_string_start(offset)
        if self.code[offset] in ')]}':
            return self._find_parens_start(offset)
        if self._is_id_char(offset):
            return self._find_word_start(offset)
        return old_offset

    def _find_primary_without_dot_start(self, offset):
        """It tries to find the undotted primary start

        It is different from `self._get_atom_start()` in that it
        follows function calls, too; such as in ``f(x)``.

        """
        last_atom = offset
        offset = self._find_last_non_space_char(last_atom)
        while offset > 0 and self.code[offset] in ')]':
            last_atom = self._find_parens_start(offset)
            offset = self._find_last_non_space_char(last_atom - 1)
        if offset >= 0 and (self.code[offset] in '"\'})]' or
                            self._is_id_char(offset)):
            atom_start = self._find_atom_start(offset)
            if not keyword.iskeyword(self.code[atom_start:offset + 1]):
                return atom_start
        return last_atom

    def _find_primary_start(self, offset):
        if offset >= len(self.code):
            offset = len(self.code) - 1
        if self.code[offset] != '.':
            offset = self._find_primary_without_dot_start(offset)
        else:
            offset = offset + 1
        while offset > 0:
            prev = self._find_last_non_space_char(offset - 1)
            if offset <= 0 or self.code[prev] != '.':
                break
            offset = self._find_primary_without_dot_start(prev - 1)
            if not self._is_id_char(offset):
                break

        return offset

    def get_primary_at(self, offset):
        offset = self._get_fixed_offset(offset)
        start, end = self.get_primary_range(offset)
        return self.raw[start:end].strip()

    def get_splitted_primary_before(self, offset):
        """returns expression, starting, starting_offset

        This function is used in `rope.codeassist.assist` function.
        """
        if offset == 0:
            return ('', '', 0)
        end = offset - 1
        word_start = self._find_atom_start(end)
        real_start = self._find_primary_start(end)
        if self.code[word_start:offset].strip() == '':
            word_start = end
        if self.code[end].isspace():
            word_start = end
        if self.code[real_start:word_start].strip() == '':
            real_start = word_start
        if real_start == word_start == end and not self._is_id_char(end):
            return ('', '', offset)
        if real_start == word_start:
            return ('', self.raw[word_start:offset], word_start)
        else:
            if self.code[end] == '.':
                return (self.raw[real_start:end], '', offset)
            last_dot_position = word_start
            if self.code[word_start] != '.':
                last_dot_position = \
                    self._find_last_non_space_char(word_start - 1)
            last_char_position = \
                self._find_last_non_space_char(last_dot_position - 1)
            if self.code[word_start].isspace():
                word_start = offset
            return (self.raw[real_start:last_char_position + 1],
                    self.raw[word_start:offset], word_start)

    def _get_line_start(self, offset):
        try:
            return self.code.rindex('\n', 0, offset + 1)
        except ValueError:
            return 0

    def _get_line_end(self, offset):
        try:
            return self.code.index('\n', offset)
        except ValueError:
            return len(self.code)

    def is_name_assigned_in_class_body(self, offset):
        word_start = self._find_word_start(offset - 1)
        word_end = self._find_word_end(offset) + 1
        if '.' in self.code[word_start:word_end]:
            return False
        line_start = self._get_line_start(word_start)
        line = self.code[line_start:word_start].strip()
        return not line and self.get_assignment_type(offset) == '='

    def is_a_class_or_function_name_in_header(self, offset):
        word_start = self._find_word_start(offset - 1)
        line_start = self._get_line_start(word_start)
        prev_word = self.code[line_start:word_start].strip()
        return prev_word in ['def', 'class']

    def _find_first_non_space_char(self, offset):
        if offset >= len(self.code):
            return len(self.code)
        while offset < len(self.code) and self.code[offset].isspace():
            if self.code[offset] == '\n':
                return offset
            offset += 1
        return offset

    def is_a_function_being_called(self, offset):
        word_end = self._find_word_end(offset) + 1
        next_char = self._find_first_non_space_char(word_end)
        return next_char < len(self.code) and \
            self.code[next_char] == '(' and \
            not self.is_a_class_or_function_name_in_header(offset)

    def _find_import_end(self, start):
        return self._get_line_end(start)

    def is_import_statement(self, offset):
        try:
            last_import = self.code.rindex('import ', 0, offset)
        except ValueError:
            return False
        return self._find_import_end(last_import + 7) >= offset

    def is_from_statement(self, offset):
        try:
            last_from = self.code.rindex('from ', 0, offset)
            from_import = self.code.index(' import ', last_from)
            from_names = from_import + 8
        except ValueError:
            return False
        from_names = self._find_first_non_space_char(from_names)
        return self._find_import_end(from_names) >= offset

    def is_from_statement_module(self, offset):
        if offset >= len(self.code) - 1:
            return False
        stmt_start = self._find_primary_start(offset)
        line_start = self._get_line_start(stmt_start)
        prev_word = self.code[line_start:stmt_start].strip()
        return prev_word == 'from'

    def is_a_name_after_from_import(self, offset):
        try:
            if len(self.code) > offset and self.code[offset] == '\n':
                line_start = self._get_line_start(offset - 1)
            else:
                line_start = self._get_line_start(offset)
            last_from = self.code.rindex('from ', line_start, offset)
            from_import = self.code.index(' import ', last_from)
            from_names = from_import + 8
        except ValueError:
            return False
        if from_names - 1 > offset:
            return False
        return self._find_import_end(from_names) >= offset

    def get_from_module(self, offset):
        try:
            last_from = self.code.rindex('from ', 0, offset)
            import_offset = self.code.index(' import ', last_from)
            end = self._find_last_non_space_char(import_offset)
            return self.get_primary_at(end)
        except ValueError:
            pass

    def is_from_aliased(self, offset):
        if not self.is_a_name_after_from_import(offset):
            return False
        try:
            end = self._find_word_end(offset)
            as_end = min(self._find_word_end(end + 1), len(self.code))
            as_start = self._find_word_start(as_end)
            if self.code[as_start:as_end + 1] == 'as':
                return True
        except ValueError:
            return False

    def get_from_aliased(self, offset):
        try:
            end = self._find_word_end(offset)
            as_ = self._find_word_end(end + 1)
            alias = self._find_word_end(as_ + 1)
            start = self._find_word_start(alias)
            return self.raw[start:alias + 1]
        except ValueError:
            pass

    def is_function_keyword_parameter(self, offset):
        word_end = self._find_word_end(offset)
        if word_end + 1 == len(self.code):
            return False
        next_char = self._find_first_non_space_char(word_end + 1)
        equals = self.code[next_char:next_char + 2]
        if equals == '==' or not equals.startswith('='):
            return False
        word_start = self._find_word_start(offset)
        prev_char = self._find_last_non_space_char(word_start - 1)
        return prev_char - 1 >= 0 and self.code[prev_char] in ',('

    def is_on_function_call_keyword(self, offset):
        stop = self._get_line_start(offset)
        if self._is_id_char(offset):
            offset = self._find_word_start(offset) - 1
        offset = self._find_last_non_space_char(offset)
        if offset <= stop or self.code[offset] not in '(,':
            return False
        parens_start = self.find_parens_start_from_inside(offset)
        return stop < parens_start

    def find_parens_start_from_inside(self, offset):
        stop = self._get_line_start(offset)
        while offset > stop:
            if self.code[offset] == '(':
                break
            if self.code[offset] != ',':
                offset = self._find_primary_start(offset)
            offset -= 1
        return max(stop, offset)

    def is_assigned_here(self, offset):
        return self.get_assignment_type(offset) is not None

    def get_assignment_type(self, offset):
        # XXX: does not handle tuple assignments
        word_end = self._find_word_end(offset)
        next_char = self._find_first_non_space_char(word_end + 1)
        single = self.code[next_char:next_char + 1]
        double = self.code[next_char:next_char + 2]
        triple = self.code[next_char:next_char + 3]
        if double not in ('==', '<=', '>=', '!='):
            for op in [single, double, triple]:
                if op.endswith('='):
                    return op

    def get_primary_range(self, offset):
        start = self._find_primary_start(offset)
        end = self._find_word_end(offset) + 1
        return (start, end)

    def get_word_range(self, offset):
        offset = max(0, offset)
        start = self._find_word_start(offset)
        end = self._find_word_end(offset) + 1
        return (start, end)

    def get_word_parens_range(self, offset, opening='(', closing=')'):
        end = self._find_word_end(offset)
        start_parens = self.code.index(opening, end)
        index = start_parens
        open_count = 0
        while index < len(self.code):
            if self.code[index] == opening:
                open_count += 1
            if self.code[index] == closing:
                open_count -= 1
            if open_count == 0:
                return (start_parens, index + 1)
            index += 1
        return (start_parens, index)

    def get_parameters(self, first, last):
        keywords = []
        args = []
        current = self._find_last_non_space_char(last - 1)
        while current > first:
            primary_start = current
            current = self._find_primary_start(current)
            while current != first and self.code[current] not in '=,':
                current = self._find_last_non_space_char(current - 1)
            primary = self.raw[current + 1:primary_start + 1].strip()
            if self.code[current] == '=':
                primary_start = current - 1
                current -= 1
                while current != first and self.code[current] not in ',':
                    current = self._find_last_non_space_char(current - 1)
                param_name = self.raw[current + 1:primary_start + 1].strip()
                keywords.append((param_name, primary))
            else:
                args.append(primary)
            current = self._find_last_non_space_char(current - 1)
        args.reverse()
        keywords.reverse()
        return args, keywords

    def is_assigned_in_a_tuple_assignment(self, offset):
        start = self._get_line_start(offset)
        end = self._get_line_end(offset)
        primary_start = self._find_primary_start(offset)
        primary_end = self._find_word_end(offset)

        prev_char_offset = self._find_last_non_space_char(primary_start - 1)
        next_char_offset = self._find_first_non_space_char(primary_end + 1)
        next_char = prev_char = ''
        if prev_char_offset >= start:
            prev_char = self.code[prev_char_offset]
        if next_char_offset < end:
            next_char = self.code[next_char_offset]
        try:
            equals_offset = self.code.index('=', start, end)
        except ValueError:
            return False
        if prev_char not in '(,' and next_char not in ',)':
            return False
        parens_start = self.find_parens_start_from_inside(offset)
        # XXX: only handling (x, y) = value
        return offset < equals_offset and \
            self.code[start:parens_start].strip() == ''

    def get_function_and_args_in_header(self, offset):
        offset = self.find_function_offset(offset)
        lparens, rparens = self.get_word_parens_range(offset)
        return self.raw[offset:rparens + 1]

    def find_function_offset(self, offset, definition='def '):
        while True:
            offset = self.code.index(definition, offset)
            if offset == 0 or not self._is_id_char(offset - 1):
                break
            offset += 1
        def_ = offset + 4
        return self._find_first_non_space_char(def_)

    def get_lambda_and_args(self, offset):
        offset = self.find_function_offset(offset, definition='lambda ')
        lparens, rparens = self.get_word_parens_range(offset, opening=' ',
                                                      closing=':')
        return self.raw[offset:rparens + 1]
