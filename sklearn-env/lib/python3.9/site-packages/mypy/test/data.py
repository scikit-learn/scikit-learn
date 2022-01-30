"""Utilities for processing .test files containing test case descriptions."""

import os.path
import os
import tempfile
import posixpath
import re
import shutil
from abc import abstractmethod
import sys

import pytest
from typing import List, Tuple, Set, Optional, Iterator, Any, Dict, NamedTuple, Union

from mypy.test.config import test_data_prefix, test_temp_dir, PREFIX

root_dir = os.path.normpath(PREFIX)

# File modify/create operation: copy module contents from source_path.
UpdateFile = NamedTuple('UpdateFile', [('module', str),
                                       ('content', str),
                                       ('target_path', str)])

# File delete operation: delete module file.
DeleteFile = NamedTuple('DeleteFile', [('module', str),
                                       ('path', str)])

FileOperation = Union[UpdateFile, DeleteFile]


def parse_test_case(case: 'DataDrivenTestCase') -> None:
    """Parse and prepare a single case from suite with test case descriptions.

    This method is part of the setup phase, just before the test case is run.
    """
    test_items = parse_test_data(case.data, case.name)
    base_path = case.suite.base_path
    if case.suite.native_sep:
        join = os.path.join
    else:
        join = posixpath.join

    out_section_missing = case.suite.required_out_section
    normalize_output = True

    files: List[Tuple[str, str]] = []  # path and contents
    output_files: List[Tuple[str, str]] = []  # path and contents for output files
    output: List[str] = []  # Regular output errors
    output2: Dict[int, List[str]] = {}  # Output errors for incremental, runs 2+
    deleted_paths: Dict[int, Set[str]] = {}  # from run number of paths
    stale_modules: Dict[int, Set[str]] = {}  # from run number to module names
    rechecked_modules: Dict[int, Set[str]] = {}  # from run number module names
    triggered: List[str] = []  # Active triggers (one line per incremental step)
    targets: Dict[int, List[str]] = {}  # Fine-grained targets (per fine-grained update)

    # Process the parsed items. Each item has a header of form [id args],
    # optionally followed by lines of text.
    item = first_item = test_items[0]
    for item in test_items[1:]:
        if item.id == 'file' or item.id == 'outfile':
            # Record an extra file needed for the test case.
            assert item.arg is not None
            contents = expand_variables('\n'.join(item.data))
            file_entry = (join(base_path, item.arg), contents)
            if item.id == 'file':
                files.append(file_entry)
            else:
                output_files.append(file_entry)
        elif item.id in ('builtins', 'builtins_py2'):
            # Use an alternative stub file for the builtins module.
            assert item.arg is not None
            mpath = join(os.path.dirname(case.file), item.arg)
            fnam = 'builtins.pyi' if item.id == 'builtins' else '__builtin__.pyi'
            with open(mpath, encoding='utf8') as f:
                files.append((join(base_path, fnam), f.read()))
        elif item.id == 'typing':
            # Use an alternative stub file for the typing module.
            assert item.arg is not None
            src_path = join(os.path.dirname(case.file), item.arg)
            with open(src_path, encoding='utf8') as f:
                files.append((join(base_path, 'typing.pyi'), f.read()))
        elif re.match(r'stale[0-9]*$', item.id):
            passnum = 1 if item.id == 'stale' else int(item.id[len('stale'):])
            assert passnum > 0
            modules = (set() if item.arg is None else {t.strip() for t in item.arg.split(',')})
            stale_modules[passnum] = modules
        elif re.match(r'rechecked[0-9]*$', item.id):
            passnum = 1 if item.id == 'rechecked' else int(item.id[len('rechecked'):])
            assert passnum > 0
            modules = (set() if item.arg is None else {t.strip() for t in item.arg.split(',')})
            rechecked_modules[passnum] = modules
        elif re.match(r'targets[0-9]*$', item.id):
            passnum = 1 if item.id == 'targets' else int(item.id[len('targets'):])
            assert passnum > 0
            reprocessed = [] if item.arg is None else [t.strip() for t in item.arg.split(',')]
            targets[passnum] = reprocessed
        elif item.id == 'delete':
            # File/directory to delete during a multi-step test case
            assert item.arg is not None
            m = re.match(r'(.*)\.([0-9]+)$', item.arg)
            assert m, 'Invalid delete section: {}'.format(item.arg)
            num = int(m.group(2))
            assert num >= 2, "Can't delete during step {}".format(num)
            full = join(base_path, m.group(1))
            deleted_paths.setdefault(num, set()).add(full)
        elif re.match(r'out[0-9]*$', item.id):
            if item.arg is None:
                args = []
            else:
                args = item.arg.split(",")

            version_check = True
            for arg in args:
                if arg == 'skip-path-normalization':
                    normalize_output = False
                if arg.startswith("version"):
                    compare_op = arg[7:9]
                    if compare_op not in {">=", "=="}:
                        raise ValueError(
                            "{}, line {}: Only >= and == version checks are currently supported"
                            .format(
                                case.file, item.line
                            )
                        )
                    version_str = arg[9:]
                    try:
                        version = tuple(int(x) for x in version_str.split("."))
                    except ValueError:
                        raise ValueError(
                            '{}, line {}: "{}" is not a valid python version'.format(
                                case.file, item.line, version_str))
                    if compare_op == ">=":
                        version_check = sys.version_info >= version
                    elif compare_op == "==":
                        if not 1 < len(version) < 4:
                            raise ValueError(
                                '{}, line {}: Only minor or patch version checks '
                                'are currently supported with "==": "{}"'.format(
                                    case.file, item.line, version_str
                                )
                            )
                        version_check = sys.version_info[:len(version)] == version
            if version_check:
                tmp_output = [expand_variables(line) for line in item.data]
                if os.path.sep == '\\' and normalize_output:
                    tmp_output = [fix_win_path(line) for line in tmp_output]
                if item.id == 'out' or item.id == 'out1':
                    output = tmp_output
                else:
                    passnum = int(item.id[len('out'):])
                    assert passnum > 1
                    output2[passnum] = tmp_output
                out_section_missing = False
        elif item.id == 'triggered' and item.arg is None:
            triggered = item.data
        else:
            raise ValueError(
                'Invalid section header {} in {} at line {}'.format(
                    item.id, case.file, item.line))

    if out_section_missing:
        raise ValueError(
            '{}, line {}: Required output section not found'.format(
                case.file, first_item.line))

    for passnum in stale_modules.keys():
        if passnum not in rechecked_modules:
            # If the set of rechecked modules isn't specified, make it the same as the set
            # of modules with a stale public interface.
            rechecked_modules[passnum] = stale_modules[passnum]
        if (passnum in stale_modules
                and passnum in rechecked_modules
                and not stale_modules[passnum].issubset(rechecked_modules[passnum])):
            raise ValueError(
                ('Stale modules after pass {} must be a subset of rechecked '
                 'modules ({}:{})').format(passnum, case.file, first_item.line))

    input = first_item.data
    expand_errors(input, output, 'main')
    for file_path, contents in files:
        expand_errors(contents.split('\n'), output, file_path)

    case.input = input
    case.output = output
    case.output2 = output2
    case.last_line = case.line + item.line + len(item.data) - 2
    case.files = files
    case.output_files = output_files
    case.expected_stale_modules = stale_modules
    case.expected_rechecked_modules = rechecked_modules
    case.deleted_paths = deleted_paths
    case.triggered = triggered or []
    case.normalize_output = normalize_output
    case.expected_fine_grained_targets = targets


class DataDrivenTestCase(pytest.Item):
    """Holds parsed data-driven test cases, and handles directory setup and teardown."""

    # Override parent member type
    parent: "DataSuiteCollector"

    input: List[str]
    output: List[str]  # Output for the first pass
    output2: Dict[int, List[str]]  # Output for runs 2+, indexed by run number

    # full path of test suite
    file = ''
    line = 0

    # (file path, file content) tuples
    files: List[Tuple[str, str]]
    expected_stale_modules: Dict[int, Set[str]]
    expected_rechecked_modules: Dict[int, Set[str]]
    expected_fine_grained_targets: Dict[int, List[str]]

    # Whether or not we should normalize the output to standardize things like
    # forward vs backward slashes in file paths for Windows vs Linux.
    normalize_output = True

    # Extra attributes used by some tests.
    last_line: int
    output_files: List[Tuple[str, str]]  # Path and contents for output files
    deleted_paths: Dict[int, Set[str]]  # Mapping run number -> paths
    triggered: List[str]  # Active triggers (one line per incremental step)

    def __init__(self,
                 parent: 'DataSuiteCollector',
                 suite: 'DataSuite',
                 file: str,
                 name: str,
                 writescache: bool,
                 only_when: str,
                 platform: Optional[str],
                 skip: bool,
                 xfail: bool,
                 data: str,
                 line: int) -> None:
        super().__init__(name, parent)
        self.suite = suite
        self.file = file
        self.writescache = writescache
        self.only_when = only_when
        if ((platform == 'windows' and sys.platform != 'win32')
                or (platform == 'posix' and sys.platform == 'win32')):
            skip = True
        self.skip = skip
        self.xfail = xfail
        self.data = data
        self.line = line
        self.old_cwd: Optional[str] = None
        self.tmpdir: Optional[tempfile.TemporaryDirectory[str]] = None

    def runtest(self) -> None:
        if self.skip:
            pytest.skip()
        # TODO: add a better error message for when someone uses skip and xfail at the same time
        elif self.xfail:
            self.add_marker(pytest.mark.xfail)
        parent = self.getparent(DataSuiteCollector)
        assert parent is not None, 'Should not happen'
        suite = parent.obj()
        suite.setup()
        try:
            suite.run_case(self)
        except Exception:
            # As a debugging aid, support copying the contents of the tmp directory somewhere
            save_dir: Optional[str] = self.config.getoption("--save-failures-to", None)
            if save_dir:
                assert self.tmpdir is not None
                target_dir = os.path.join(save_dir, os.path.basename(self.tmpdir.name))
                print("Copying data from test {} to {}".format(self.name, target_dir))
                if not os.path.isabs(target_dir):
                    assert self.old_cwd
                    target_dir = os.path.join(self.old_cwd, target_dir)
                shutil.copytree(self.tmpdir.name, target_dir)
            raise

    def setup(self) -> None:
        parse_test_case(case=self)
        self.old_cwd = os.getcwd()
        self.tmpdir = tempfile.TemporaryDirectory(prefix='mypy-test-')
        os.chdir(self.tmpdir.name)
        os.mkdir(test_temp_dir)

        # Precalculate steps for find_steps()
        steps: Dict[int, List[FileOperation]] = {}

        for path, content in self.files:
            m = re.match(r'.*\.([0-9]+)$', path)
            if m:
                # Skip writing subsequent incremental steps - rather
                # store them as operations.
                num = int(m.group(1))
                assert num >= 2
                target_path = re.sub(r'\.[0-9]+$', '', path)
                module = module_from_path(target_path)
                operation = UpdateFile(module, content, target_path)
                steps.setdefault(num, []).append(operation)
            else:
                # Write the first incremental steps
                dir = os.path.dirname(path)
                os.makedirs(dir, exist_ok=True)
                with open(path, 'w', encoding='utf8') as f:
                    f.write(content)

        for num, paths in self.deleted_paths.items():
            assert num >= 2
            for path in paths:
                module = module_from_path(path)
                steps.setdefault(num, []).append(DeleteFile(module, path))
        max_step = max(steps) if steps else 2
        self.steps = [steps.get(num, []) for num in range(2, max_step + 1)]

    def teardown(self) -> None:
        assert self.old_cwd is not None and self.tmpdir is not None, \
            "test was not properly set up"
        os.chdir(self.old_cwd)
        try:
            self.tmpdir.cleanup()
        except OSError:
            pass
        self.old_cwd = None
        self.tmpdir = None

    def reportinfo(self) -> Tuple[str, int, str]:
        return self.file, self.line, self.name

    def repr_failure(self, excinfo: Any, style: Optional[Any] = None) -> str:
        if excinfo.errisinstance(SystemExit):
            # We assume that before doing exit() (which raises SystemExit) we've printed
            # enough context about what happened so that a stack trace is not useful.
            # In particular, uncaught exceptions during semantic analysis or type checking
            # call exit() and they already print out a stack trace.
            excrepr = excinfo.exconly()
        else:
            self.parent._prunetraceback(excinfo)
            excrepr = excinfo.getrepr(style='short')

        return "data: {}:{}:\n{}".format(self.file, self.line, excrepr)

    def find_steps(self) -> List[List[FileOperation]]:
        """Return a list of descriptions of file operations for each incremental step.

        The first list item corresponds to the first incremental step, the second for the
        second step, etc. Each operation can either be a file modification/creation (UpdateFile)
        or deletion (DeleteFile).

        Defaults to having two steps if there aern't any operations.
        """
        return self.steps


def module_from_path(path: str) -> str:
    path = re.sub(r'\.pyi?$', '', path)
    # We can have a mix of Unix-style and Windows-style separators.
    parts = re.split(r'[/\\]', path)
    del parts[0]
    module = '.'.join(parts)
    module = re.sub(r'\.__init__$', '', module)
    return module


class TestItem:
    """Parsed test caseitem.

    An item is of the form
      [id arg]
      .. data ..
    """

    id = ""
    arg: Optional[str] = ""

    # Text data, array of 8-bit strings
    data: List[str]

    file = ''
    line = 0  # Line number in file

    def __init__(self, id: str, arg: Optional[str], data: List[str],
                 line: int) -> None:
        self.id = id
        self.arg = arg
        self.data = data
        self.line = line


def parse_test_data(raw_data: str, name: str) -> List[TestItem]:
    """Parse a list of lines that represent a sequence of test items."""

    lines = ['', '[case ' + name + ']'] + raw_data.split('\n')
    ret: List[TestItem] = []
    data: List[str] = []

    id: Optional[str] = None
    arg: Optional[str] = None

    i = 0
    i0 = 0
    while i < len(lines):
        s = lines[i].strip()

        if lines[i].startswith('[') and s.endswith(']'):
            if id:
                data = collapse_line_continuation(data)
                data = strip_list(data)
                ret.append(TestItem(id, arg, strip_list(data), i0 + 1))

            i0 = i
            id = s[1:-1]
            arg = None
            if ' ' in id:
                arg = id[id.index(' ') + 1:]
                id = id[:id.index(' ')]
            data = []
        elif lines[i].startswith('\\['):
            data.append(lines[i][1:])
        elif not lines[i].startswith('--'):
            data.append(lines[i])
        elif lines[i].startswith('----'):
            data.append(lines[i][2:])
        i += 1

    # Process the last item.
    if id:
        data = collapse_line_continuation(data)
        data = strip_list(data)
        ret.append(TestItem(id, arg, data, i0 + 1))

    return ret


def strip_list(l: List[str]) -> List[str]:
    """Return a stripped copy of l.

    Strip whitespace at the end of all lines, and strip all empty
    lines from the end of the array.
    """

    r: List[str] = []
    for s in l:
        # Strip spaces at end of line
        r.append(re.sub(r'\s+$', '', s))

    while len(r) > 0 and r[-1] == '':
        r.pop()

    return r


def collapse_line_continuation(l: List[str]) -> List[str]:
    r: List[str] = []
    cont = False
    for s in l:
        ss = re.sub(r'\\$', '', s)
        if cont:
            r[-1] += re.sub('^ +', '', ss)
        else:
            r.append(ss)
        cont = s.endswith('\\')
    return r


def expand_variables(s: str) -> str:
    return s.replace('<ROOT>', root_dir)


def expand_errors(input: List[str], output: List[str], fnam: str) -> None:
    """Transform comments such as '# E: message' or
    '# E:3: message' in input.

    The result is lines like 'fnam:line: error: message'.
    """

    for i in range(len(input)):
        # The first in the split things isn't a comment
        for possible_err_comment in input[i].split(' # ')[1:]:
            m = re.search(
                r'^([ENW]):((?P<col>\d+):)? (?P<message>.*)$',
                possible_err_comment.strip())
            if m:
                if m.group(1) == 'E':
                    severity = 'error'
                elif m.group(1) == 'N':
                    severity = 'note'
                elif m.group(1) == 'W':
                    severity = 'warning'
                col = m.group('col')
                message = m.group('message')
                message = message.replace('\\#', '#')  # adds back escaped # character
                if col is None:
                    output.append(
                        '{}:{}: {}: {}'.format(fnam, i + 1, severity, message))
                else:
                    output.append('{}:{}:{}: {}: {}'.format(
                        fnam, i + 1, col, severity, message))


def fix_win_path(line: str) -> str:
    r"""Changes Windows paths to Linux paths in error messages.

    E.g. foo\bar.py -> foo/bar.py.
    """
    line = line.replace(root_dir, root_dir.replace('\\', '/'))
    m = re.match(r'^([\S/]+):(\d+:)?(\s+.*)', line)
    if not m:
        return line
    else:
        filename, lineno, message = m.groups()
        return '{}:{}{}'.format(filename.replace('\\', '/'),
                                lineno or '', message)


def fix_cobertura_filename(line: str) -> str:
    r"""Changes filename paths to Linux paths in Cobertura output files.

    E.g. filename="pkg\subpkg\a.py" -> filename="pkg/subpkg/a.py".
    """
    m = re.search(r'<class .* filename="(?P<filename>.*?)"', line)
    if not m:
        return line
    return '{}{}{}'.format(line[:m.start(1)],
                           m.group('filename').replace('\\', '/'),
                           line[m.end(1):])


##
#
# pytest setup
#
##


# This function name is special to pytest.  See
# https://docs.pytest.org/en/latest/reference.html#initialization-hooks
def pytest_addoption(parser: Any) -> None:
    group = parser.getgroup('mypy')
    group.addoption('--update-data', action='store_true', default=False,
                    help='Update test data to reflect actual output'
                         ' (supported only for certain tests)')
    group.addoption('--save-failures-to', default=None,
                    help='Copy the temp directories from failing tests to a target directory')
    group.addoption('--mypy-verbose', action='count',
                    help='Set the verbose flag when creating mypy Options')
    group.addoption('--mypyc-showc', action='store_true', default=False,
                    help='Display C code on mypyc test failures')


# This function name is special to pytest.  See
# http://doc.pytest.org/en/latest/writing_plugins.html#collection-hooks
def pytest_pycollect_makeitem(collector: Any, name: str,
                              obj: object) -> 'Optional[Any]':
    """Called by pytest on each object in modules configured in conftest.py files.

    collector is pytest.Collector, returns Optional[pytest.Class]
    """
    if isinstance(obj, type):
        # Only classes derived from DataSuite contain test cases, not the DataSuite class itself
        if issubclass(obj, DataSuite) and obj is not DataSuite:
            # Non-None result means this obj is a test case.
            # The collect method of the returned DataSuiteCollector instance will be called later,
            # with self.obj being obj.
            return DataSuiteCollector.from_parent(  # type: ignore[no-untyped-call]
                parent=collector, name=name,
            )
    return None


def split_test_cases(parent: 'DataFileCollector', suite: 'DataSuite',
                     file: str) -> Iterator['DataDrivenTestCase']:
    """Iterate over raw test cases in file, at collection time, ignoring sub items.

    The collection phase is slow, so any heavy processing should be deferred to after
    uninteresting tests are filtered (when using -k PATTERN switch).
    """
    with open(file, encoding='utf-8') as f:
        data = f.read()
    # number of groups in the below regex
    NUM_GROUPS = 7
    cases = re.split(r'^\[case ([a-zA-Z_0-9]+)'
                     r'(-writescache)?'
                     r'(-only_when_cache|-only_when_nocache)?'
                     r'(-posix|-windows)?'
                     r'(-skip)?'
                     r'(-xfail)?'
                     r'\][ \t]*$\n',
                     data,
                     flags=re.DOTALL | re.MULTILINE)
    line_no = cases[0].count('\n') + 1
    test_names = set()
    for i in range(1, len(cases), NUM_GROUPS):
        name, writescache, only_when, platform_flag, skip, xfail, data = cases[i:i + NUM_GROUPS]
        if name in test_names:
            raise RuntimeError('Found a duplicate test name "{}" in {} on line {}'.format(
                name, parent.name, line_no,
            ))
        platform = platform_flag[1:] if platform_flag else None
        yield DataDrivenTestCase.from_parent(
            parent=parent,
            suite=suite,
            file=file,
            name=add_test_name_suffix(name, suite.test_name_suffix),
            writescache=bool(writescache),
            only_when=only_when,
            platform=platform,
            skip=bool(skip),
            xfail=bool(xfail),
            data=data,
            line=line_no,
        )
        line_no += data.count('\n') + 1

        # Record existing tests to prevent duplicates:
        test_names.update({name})


class DataSuiteCollector(pytest.Class):
    def collect(self) -> Iterator['DataFileCollector']:
        """Called by pytest on each of the object returned from pytest_pycollect_makeitem"""

        # obj is the object for which pytest_pycollect_makeitem returned self.
        suite: DataSuite = self.obj

        assert os.path.isdir(suite.data_prefix), \
            'Test data prefix ({}) not set correctly'.format(suite.data_prefix)

        for data_file in suite.files:
            yield DataFileCollector.from_parent(parent=self, name=data_file)


class DataFileCollector(pytest.Collector):
    """Represents a single `.test` data driven test file.

    More context: https://github.com/python/mypy/issues/11662
    """
    parent: DataSuiteCollector

    @classmethod  # We have to fight with pytest here:
    def from_parent(  # type: ignore[override]
        cls,
        parent: DataSuiteCollector,
        *,
        name: str,
    ) -> 'DataFileCollector':
        return super().from_parent(parent, name=name)

    def collect(self) -> Iterator['DataDrivenTestCase']:
        yield from split_test_cases(
            parent=self,
            suite=self.parent.obj,
            file=os.path.join(self.parent.obj.data_prefix, self.name),
        )


def add_test_name_suffix(name: str, suffix: str) -> str:
    # Find magic suffix of form "-foobar" (used for things like "-skip").
    m = re.search(r'-[-A-Za-z0-9]+$', name)
    if m:
        # Insert suite-specific test name suffix before the magic suffix
        # which must be the last thing in the test case name since we
        # are using endswith() checks.
        magic_suffix = m.group(0)
        return name[:-len(magic_suffix)] + suffix + magic_suffix
    else:
        return name + suffix


def is_incremental(testcase: DataDrivenTestCase) -> bool:
    return 'incremental' in testcase.name.lower() or 'incremental' in testcase.file


def has_stable_flags(testcase: DataDrivenTestCase) -> bool:
    if any(re.match(r'# flags[2-9]:', line) for line in testcase.input):
        return False
    for filename, contents in testcase.files:
        if os.path.basename(filename).startswith('mypy.ini.'):
            return False
    return True


class DataSuite:
    # option fields - class variables
    files: List[str]

    base_path = test_temp_dir

    # Allow external users of the test code to override the data prefix
    data_prefix = test_data_prefix

    required_out_section = False

    native_sep = False

    # Name suffix automatically added to each test case in the suite (can be
    # used to distinguish test cases in suites that share data files)
    test_name_suffix = ''

    def setup(self) -> None:
        """Setup fixtures (ad-hoc)"""
        pass

    @abstractmethod
    def run_case(self, testcase: DataDrivenTestCase) -> None:
        raise NotImplementedError
