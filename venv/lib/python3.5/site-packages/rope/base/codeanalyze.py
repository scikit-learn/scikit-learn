import bisect
import re
import token
import tokenize


class ChangeCollector(object):

    def __init__(self, text):
        self.text = text
        self.changes = []

    def add_change(self, start, end, new_text=None):
        if new_text is None:
            new_text = self.text[start:end]
        self.changes.append((start, end, new_text))

    def get_changed(self):
        if not self.changes:
            return None

        self.changes.sort(key=lambda x: x[:2])
        pieces = []
        last_changed = 0
        for change in self.changes:
            start, end, text = change
            pieces.append(self.text[last_changed:start] + text)
            last_changed = end
        if last_changed < len(self.text):
            pieces.append(self.text[last_changed:])
        result = ''.join(pieces)
        if result != self.text:
            return result


class SourceLinesAdapter(object):
    """Adapts source to Lines interface

    Note: The creation of this class is expensive.
    """

    def __init__(self, source_code):
        self.code = source_code
        self.starts = None
        self._initialize_line_starts()

    def _initialize_line_starts(self):
        self.starts = []
        self.starts.append(0)
        try:
            i = 0
            while True:
                i = self.code.index('\n', i) + 1
                self.starts.append(i)
        except ValueError:
            pass
        self.starts.append(len(self.code) + 1)

    def get_line(self, lineno):
        return self.code[self.starts[lineno - 1]:
                         self.starts[lineno] - 1]

    def length(self):
        return len(self.starts) - 1

    def get_line_number(self, offset):
        return bisect.bisect(self.starts, offset)

    def get_line_start(self, lineno):
        return self.starts[lineno - 1]

    def get_line_end(self, lineno):
        return self.starts[lineno] - 1


class ArrayLinesAdapter(object):

    def __init__(self, lines):
        self.lines = lines

    def get_line(self, line_number):
        return self.lines[line_number - 1]

    def length(self):
        return len(self.lines)


class LinesToReadline(object):

    def __init__(self, lines, start):
        self.lines = lines
        self.current = start

    def readline(self):
        if self.current <= self.lines.length():
            self.current += 1
            return self.lines.get_line(self.current - 1) + '\n'
        return ''

    def __call__(self):
        return self.readline()


class _CustomGenerator(object):

    def __init__(self, lines):
        self.lines = lines
        self.in_string = ''
        self.open_count = 0
        self.continuation = False

    def __call__(self):
        size = self.lines.length()
        result = []
        i = 1
        while i <= size:
            while i <= size and not self.lines.get_line(i).strip():
                i += 1
            if i <= size:
                start = i
                while True:
                    line = self.lines.get_line(i)
                    self._analyze_line(line)
                    if not (self.continuation or self.open_count or
                            self.in_string) or i == size:
                        break
                    i += 1
                result.append((start, i))
                i += 1
        return result

    # Matches all backslashes before the token, to detect escaped quotes
    _main_tokens = re.compile(r'(\\*)((\'\'\'|"""|\'|")|#|\[|\]|\{|\}|\(|\))')

    def _analyze_line(self, line):
        token = None
        for match in self._main_tokens.finditer(line):
            prefix = match.group(1)
            token = match.group(2)
            # Skip any tokens which are escaped
            if len(prefix) % 2 == 1:
                continue
            if token in ["'''", '"""', "'", '"']:
                if not self.in_string:
                    self.in_string = token
                elif self.in_string == token:
                    self.in_string = ''
            if self.in_string:
                continue
            if token == '#':
                break
            if token in '([{':
                self.open_count += 1
            elif token in ')]}':
                self.open_count -= 1
        if line and token != '#' and line.endswith('\\'):
            self.continuation = True
        else:
            self.continuation = False


def custom_generator(lines):
    return _CustomGenerator(lines)()


class LogicalLineFinder(object):

    def __init__(self, lines):
        self.lines = lines

    def logical_line_in(self, line_number):
        indents = count_line_indents(self.lines.get_line(line_number))
        tries = 0
        while True:
            block_start = get_block_start(self.lines, line_number, indents)
            try:
                return self._block_logical_line(block_start, line_number)
            except IndentationError as e:
                tries += 1
                if tries == 5:
                    raise e
                lineno = e.lineno + block_start - 1
                indents = count_line_indents(self.lines.get_line(lineno))

    def generate_starts(self, start_line=1, end_line=None):
        for start, end in self.generate_regions(start_line, end_line):
            yield start

    def generate_regions(self, start_line=1, end_line=None):
        # XXX: `block_start` should be at a better position!
        block_start = 1
        readline = LinesToReadline(self.lines, block_start)
        try:
            for start, end in self._logical_lines(readline):
                real_start = start + block_start - 1
                real_start = self._first_non_blank(real_start)
                if end_line is not None and real_start >= end_line:
                    break
                real_end = end + block_start - 1
                if real_start >= start_line:
                    yield (real_start, real_end)
        except tokenize.TokenError:
            pass

    def _block_logical_line(self, block_start, line_number):
        readline = LinesToReadline(self.lines, block_start)
        shifted = line_number - block_start + 1
        region = self._calculate_logical(readline, shifted)
        start = self._first_non_blank(region[0] + block_start - 1)
        if region[1] is None:
            end = self.lines.length()
        else:
            end = region[1] + block_start - 1
        return start, end

    def _calculate_logical(self, readline, line_number):
        last_end = 1
        try:
            for start, end in self._logical_lines(readline):
                if line_number <= end:
                    return (start, end)
                last_end = end + 1
        except tokenize.TokenError as e:
            current = e.args[1][0]
            return (last_end, max(last_end, current - 1))
        return (last_end, None)

    def _logical_lines(self, readline):
        last_end = 1
        for current_token in tokenize.generate_tokens(readline):
            current = current_token[2][0]
            if current_token[0] == token.NEWLINE:
                yield (last_end, current)
                last_end = current + 1

    def _first_non_blank(self, line_number):
        current = line_number
        while current < self.lines.length():
            line = self.lines.get_line(current).strip()
            if line and not line.startswith('#'):
                return current
            current += 1
        return current


def tokenizer_generator(lines):
    return LogicalLineFinder(lines).generate_regions()


class CachingLogicalLineFinder(object):

    def __init__(self, lines, generate=custom_generator):
        self.lines = lines
        self._generate = generate

    _starts = None

    @property
    def starts(self):
        if self._starts is None:
            self._init_logicals()
        return self._starts

    _ends = None

    @property
    def ends(self):
        if self._ends is None:
            self._init_logicals()
        return self._ends

    def _init_logicals(self):
        """Should initialize _starts and _ends attributes"""
        size = self.lines.length() + 1
        self._starts = [None] * size
        self._ends = [None] * size
        for start, end in self._generate(self.lines):
            self._starts[start] = True
            self._ends[end] = True

    def logical_line_in(self, line_number):
        start = line_number
        while start > 0 and not self.starts[start]:
            start -= 1
        if start == 0:
            try:
                start = self.starts.index(True, line_number)
            except ValueError:
                return (line_number, line_number)
        return (start, self.ends.index(True, start))

    def generate_starts(self, start_line=1, end_line=None):
        if end_line is None:
            end_line = self.lines.length()
        for index in range(start_line, end_line):
            if self.starts[index]:
                yield index


def get_block_start(lines, lineno, maximum_indents=80):
    """Approximate block start"""
    pattern = get_block_start_patterns()
    for i in range(lineno, 0, -1):
        match = pattern.search(lines.get_line(i))
        if match is not None and \
           count_line_indents(lines.get_line(i)) <= maximum_indents:
            striped = match.string.lstrip()
            # Maybe we're in a list comprehension or generator expression
            if i > 1 and striped.startswith('if') or striped.startswith('for'):
                bracs = 0
                for j in range(i, min(i + 5, lines.length() + 1)):
                    for c in lines.get_line(j):
                        if c == '#':
                            break
                        if c in '[(':
                            bracs += 1
                        if c in ')]':
                            bracs -= 1
                            if bracs < 0:
                                break
                    if bracs < 0:
                        break
                if bracs < 0:
                    continue
            return i
    return 1


_block_start_pattern = None


def get_block_start_patterns():
    global _block_start_pattern
    if not _block_start_pattern:
        pattern = '^\\s*(((def|class|if|elif|except|for|while|with)\\s)|'\
                  '((try|else|finally|except)\\s*:))'
        _block_start_pattern = re.compile(pattern, re.M)
    return _block_start_pattern


def count_line_indents(line):
    indents = 0
    for char in line:
        if char == ' ':
            indents += 1
        elif char == '\t':
            indents += 8
        else:
            return indents
    return 0


def get_string_pattern():
    start = r'(\b[uU]?[rR]?)?'
    longstr = r'%s"""(\\.|"(?!"")|\\\n|[^"\\])*"""' % start
    shortstr = r'%s"(\\.|\\\n|[^"\\])*"' % start
    return '|'.join([longstr, longstr.replace('"', "'"),
                     shortstr, shortstr.replace('"', "'")])


def get_comment_pattern():
    return r'#[^\n]*'
