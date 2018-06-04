"""Project file system commands.

This modules implements file system operations used by rope.  Different
version control systems can be supported by implementing the interface
provided by `FileSystemCommands` class.  See `SubversionCommands` and
`MercurialCommands` for example.

"""
import os
import shutil
import subprocess

import rope.base.utils.pycompat as pycompat

try:
    unicode
except NameError:
    unicode = str

def create_fscommands(root):
    dirlist = os.listdir(root)
    commands = {'.hg': MercurialCommands,
                '.svn': SubversionCommands,
                '.git': GITCommands,
                '_svn': SubversionCommands,
                '_darcs': DarcsCommands}
    for key in commands:
        if key in dirlist:
            try:
                return commands[key](root)
            except (ImportError, OSError):
                pass
    return FileSystemCommands()


class FileSystemCommands(object):

    def create_file(self, path):
        open(path, 'w').close()

    def create_folder(self, path):
        os.mkdir(path)

    def move(self, path, new_location):
        shutil.move(path, new_location)

    def remove(self, path):
        if os.path.isfile(path):
            os.remove(path)
        else:
            shutil.rmtree(path)

    def write(self, path, data):
        file_ = open(path, 'wb')
        try:
            file_.write(data)
        finally:
            file_.close()


class SubversionCommands(object):

    def __init__(self, *args):
        self.normal_actions = FileSystemCommands()
        import pysvn
        self.client = pysvn.Client()

    def create_file(self, path):
        self.normal_actions.create_file(path)
        self.client.add(path, force=True)

    def create_folder(self, path):
        self.normal_actions.create_folder(path)
        self.client.add(path, force=True)

    def move(self, path, new_location):
        self.client.move(path, new_location, force=True)

    def remove(self, path):
        self.client.remove(path, force=True)

    def write(self, path, data):
        self.normal_actions.write(path, data)


class MercurialCommands(object):

    def __init__(self, root):
        self.hg = self._import_mercurial()
        self.normal_actions = FileSystemCommands()
        try:
            self.ui = self.hg.ui.ui(
                verbose=False, debug=False, quiet=True,
                interactive=False, traceback=False, report_untrusted=False)
        except:
            self.ui = self.hg.ui.ui()
            self.ui.setconfig('ui', 'interactive', 'no')
            self.ui.setconfig('ui', 'debug', 'no')
            self.ui.setconfig('ui', 'traceback', 'no')
            self.ui.setconfig('ui', 'verbose', 'no')
            self.ui.setconfig('ui', 'report_untrusted', 'no')
            self.ui.setconfig('ui', 'quiet', 'yes')

        self.repo = self.hg.hg.repository(self.ui, root)

    def _import_mercurial(self):
        import mercurial.commands
        import mercurial.hg
        import mercurial.ui
        return mercurial

    def create_file(self, path):
        self.normal_actions.create_file(path)
        self.hg.commands.add(self.ui, self.repo, path)

    def create_folder(self, path):
        self.normal_actions.create_folder(path)

    def move(self, path, new_location):
        self.hg.commands.rename(self.ui, self.repo, path,
                                new_location, after=False)

    def remove(self, path):
        self.hg.commands.remove(self.ui, self.repo, path)

    def write(self, path, data):
        self.normal_actions.write(path, data)


class GITCommands(object):

    def __init__(self, root):
        self.root = root
        self._do(['version'])
        self.normal_actions = FileSystemCommands()

    def create_file(self, path):
        self.normal_actions.create_file(path)
        self._do(['add', self._in_dir(path)])

    def create_folder(self, path):
        self.normal_actions.create_folder(path)

    def move(self, path, new_location):
        self._do(['mv', self._in_dir(path), self._in_dir(new_location)])

    def remove(self, path):
        self._do(['rm', self._in_dir(path)])

    def write(self, path, data):
        # XXX: should we use ``git add``?
        self.normal_actions.write(path, data)

    def _do(self, args):
        _execute(['git'] + args, cwd=self.root)

    def _in_dir(self, path):
        if path.startswith(self.root):
            return path[len(self.root) + 1:]
        return self.root


class DarcsCommands(object):

    def __init__(self, root):
        self.root = root
        self.normal_actions = FileSystemCommands()

    def create_file(self, path):
        self.normal_actions.create_file(path)
        self._do(['add', path])

    def create_folder(self, path):
        self.normal_actions.create_folder(path)
        self._do(['add', path])

    def move(self, path, new_location):
        self._do(['mv', path, new_location])

    def remove(self, path):
        self.normal_actions.remove(path)

    def write(self, path, data):
        self.normal_actions.write(path, data)

    def _do(self, args):
        _execute(['darcs'] + args, cwd=self.root)


def _execute(args, cwd=None):
    process = subprocess.Popen(args, cwd=cwd, stdout=subprocess.PIPE)
    process.wait()
    return process.returncode


def unicode_to_file_data(contents, encoding=None):
    if not isinstance(contents, unicode):
        return contents
    if encoding is None:
        encoding = read_str_coding(contents)
    if encoding is not None:
        return contents.encode(encoding)
    try:
        return contents.encode()
    except UnicodeEncodeError:
        return contents.encode('utf-8')


def file_data_to_unicode(data, encoding=None):
    result = _decode_data(data, encoding)
    if '\r' in result:
        result = result.replace('\r\n', '\n').replace('\r', '\n')
    return result


def _decode_data(data, encoding):
    if isinstance(data, unicode):
        return data
    if encoding is None:
        encoding = read_str_coding(data)
    if encoding is None:
        # there is no encoding tip, we need to guess.
        # PEP263 says that "encoding not explicitly defined" means it is ascii,
        # but we will use utf8 instead since utf8 fully covers ascii and btw is
        # the only non-latin sane encoding.
        encoding = 'utf-8'
    try:
        return data.decode(encoding)
    except (UnicodeError, LookupError):
        # fallback to latin1: it should never fail
        return data.decode('latin1')


def read_file_coding(path):
    file = open(path, 'b')
    count = 0
    result = []
    while True:
        current = file.read(10)
        if not current:
            break
        count += current.count('\n')
        result.append(current)
    file.close()
    return _find_coding(''.join(result))


def read_str_coding(source):
    if type(source) == bytes:
        newline = b'\n'
    else:
        newline = '\n'
    #try:
    #    source = source.decode("utf-8")
    #except AttributeError:
    #    pass
    try:
        first = source.index(newline) + 1
        second = source.index(newline, first) + 1
    except ValueError:
        second = len(source)
    return _find_coding(source[:second])


def _find_coding(text):
    if isinstance(text, pycompat.str):
        text = text.encode('utf-8')
    coding = b'coding'
    to_chr = chr if pycompat.PY3 else lambda x: x
    try:
        start = text.index(coding) + len(coding)
        if text[start] not in b'=:':
            return
        start += 1
        while start < len(text) and to_chr(text[start]).isspace():
            start += 1
        end = start
        while end < len(text):
            c = text[end]
            if not to_chr(c).isalnum() and c not in b'-_':
                break
            end += 1
        result = text[start:end]
        if isinstance(result, bytes):
            result = result.decode('utf-8')
        return result
    except ValueError:
        pass
