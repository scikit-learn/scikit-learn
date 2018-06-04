import os
import shutil
import sys
import warnings

import rope.base.fscommands
import rope.base.resourceobserver as resourceobserver
import rope.base.utils.pycompat as pycompat
from rope.base import exceptions, taskhandle, prefs, history, pycore, utils
from rope.base.exceptions import ModuleNotFoundError
from rope.base.resources import File, Folder, _ResourceMatcher

try:
    import pickle
except ImportError:
    import cPickle as pickle


class _Project(object):

    def __init__(self, fscommands):
        self.observers = []
        self.fscommands = fscommands
        self.prefs = prefs.Prefs()
        self.data_files = _DataFiles(self)
        self._custom_source_folders = []

    def get_resource(self, resource_name):
        """Get a resource in a project.

        `resource_name` is the path of a resource in a project.  It is
        the path of a resource relative to project root.  Project root
        folder address is an empty string.  If the resource does not
        exist a `exceptions.ResourceNotFound` exception would be
        raised.  Use `get_file()` and `get_folder()` when you need to
        get nonexistent `Resource`\s.

        """
        path = self._get_resource_path(resource_name)
        if not os.path.exists(path):
            raise exceptions.ResourceNotFoundError(
                'Resource <%s> does not exist' % resource_name)
        elif os.path.isfile(path):
            return File(self, resource_name)
        elif os.path.isdir(path):
            return Folder(self, resource_name)
        else:
            raise exceptions.ResourceNotFoundError('Unknown resource '
                                                   + resource_name)

    def get_module(self, name, folder=None):
        """Returns a `PyObject` if the module was found."""
        # check if this is a builtin module
        pymod = self.pycore.builtin_module(name)
        if pymod is not None:
            return pymod
        module = self.find_module(name, folder)
        if module is None:
            raise ModuleNotFoundError('Module %s not found' % name)
        return self.pycore.resource_to_pyobject(module)

    def get_python_path_folders(self):
        result = []
        for src in self.prefs.get('python_path', []) + sys.path:
            try:
                src_folder = get_no_project().get_resource(src)
                result.append(src_folder)
            except exceptions.ResourceNotFoundError:
                pass
        return result

    # INFO: It was decided not to cache source folders, since:
    #  - Does not take much time when the root folder contains
    #    packages, that is most of the time
    #  - We need a separate resource observer; `self.observer`
    #    does not get notified about module and folder creations
    def get_source_folders(self):
        """Returns project source folders"""
        if self.root is None:
            return []
        result = list(self._custom_source_folders)
        result.extend(self.pycore._find_source_folders(self.root))
        return result

    def validate(self, folder):
        """Validate files and folders contained in this folder

        It validates all of the files and folders contained in this
        folder if some observers are interested in them.

        """
        for observer in list(self.observers):
            observer.validate(folder)

    def add_observer(self, observer):
        """Register a `ResourceObserver`

        See `FilteredResourceObserver`.
        """
        self.observers.append(observer)

    def remove_observer(self, observer):
        """Remove a registered `ResourceObserver`"""
        if observer in self.observers:
            self.observers.remove(observer)

    def do(self, changes, task_handle=taskhandle.NullTaskHandle()):
        """Apply the changes in a `ChangeSet`

        Most of the time you call this function for committing the
        changes for a refactoring.
        """
        self.history.do(changes, task_handle=task_handle)

    def get_pymodule(self, resource, force_errors=False):
        return self.pycore.resource_to_pyobject(resource, force_errors)

    def get_pycore(self):
        return self.pycore

    def get_file(self, path):
        """Get the file with `path` (it may not exist)"""
        return File(self, path)

    def get_folder(self, path):
        """Get the folder with `path` (it may not exist)"""
        return Folder(self, path)

    def get_prefs(self):
        return self.prefs

    def get_relative_module(self, name, folder, level):
        module = self.find_relative_module(name, folder, level)
        if module is None:
            raise ModuleNotFoundError('Module %s not found' % name)
        return self.pycore.resource_to_pyobject(module)

    def find_module(self, modname, folder=None):
        """Returns a resource corresponding to the given module

        returns None if it can not be found
        """
        for src in self.get_source_folders():
            module = _find_module_in_folder(src, modname)
            if module is not None:
                return module
        for src in self.get_python_path_folders():
            module = _find_module_in_folder(src, modname)
            if module is not None:
                return module
        if folder is not None:
            module = _find_module_in_folder(folder, modname)
            if module is not None:
                return module
        return None

    def find_relative_module(self, modname, folder, level):
        for i in range(level - 1):
            folder = folder.parent
        if modname == '':
            return folder
        else:
            return _find_module_in_folder(folder, modname)

    def is_ignored(self, resource):
        return False

    def _get_resource_path(self, name):
        pass

    @property
    @utils.saveit
    def history(self):
        return history.History(self)

    @property
    @utils.saveit
    def pycore(self):
        return pycore.PyCore(self)

    def close(self):
        warnings.warn('Cannot close a NoProject',
                      DeprecationWarning, stacklevel=2)

    ropefolder = None


class Project(_Project):
    """A Project containing files and folders"""

    def __init__(self, projectroot, fscommands=None,
                 ropefolder='.ropeproject', **prefs):
        """A rope project

        :parameters:
            - `projectroot`: The address of the root folder of the project
            - `fscommands`: Implements the file system operations used
              by rope; have a look at `rope.base.fscommands`
            - `ropefolder`: The name of the folder in which rope stores
              project configurations and data.  Pass `None` for not using
              such a folder at all.
            - `prefs`: Specify project preferences.  These values
              overwrite config file preferences.

        """
        if projectroot != '/':
            projectroot = _realpath(projectroot).rstrip('/\\')
        self._address = projectroot
        self._ropefolder_name = ropefolder
        if not os.path.exists(self._address):
            os.mkdir(self._address)
        elif not os.path.isdir(self._address):
            raise exceptions.RopeError('Project root exists and'
                                       ' is not a directory')
        if fscommands is None:
            fscommands = rope.base.fscommands.create_fscommands(self._address)
        super(Project, self).__init__(fscommands)
        self.ignored = _ResourceMatcher()
        self.file_list = _FileListCacher(self)
        self.prefs.add_callback('ignored_resources', self.ignored.set_patterns)
        if ropefolder is not None:
            self.prefs['ignored_resources'] = [ropefolder]
        self._init_prefs(prefs)
        self._init_source_folders()

    @utils.deprecated('Delete once deprecated functions are gone')
    def _init_source_folders(self):
        for path in self.prefs.get('source_folders', []):
            folder = self.get_resource(path)
            self._custom_source_folders.append(folder)

    def get_files(self):
        return self.file_list.get_files()

    def get_python_files(self):
        """Returns all python files available in the project"""
        return [resource for resource in self.get_files()
                if self.pycore.is_python_file(resource)]

    def _get_resource_path(self, name):
        return os.path.join(self._address, *name.split('/'))

    def _init_ropefolder(self):
        if self.ropefolder is not None:
            if not self.ropefolder.exists():
                self._create_recursively(self.ropefolder)
            if not self.ropefolder.has_child('config.py'):
                config = self.ropefolder.create_file('config.py')
                config.write(self._default_config())

    def _create_recursively(self, folder):
        if folder.parent != self.root and not folder.parent.exists():
            self._create_recursively(folder.parent)
        folder.create()

    def _init_prefs(self, prefs):
        run_globals = {}
        if self.ropefolder is not None:
            config = self.get_file(self.ropefolder.path + '/config.py')
            run_globals.update({'__name__': '__main__',
                                '__builtins__': __builtins__,
                                '__file__': config.real_path})
            if config.exists():
                config = self.ropefolder.get_child('config.py')
                pycompat.execfile(config.real_path, run_globals)
            else:
                exec(self._default_config(), run_globals)
            if 'set_prefs' in run_globals:
                run_globals['set_prefs'](self.prefs)
        for key, value in prefs.items():
            self.prefs[key] = value
        self._init_other_parts()
        self._init_ropefolder()
        if 'project_opened' in run_globals:
            run_globals['project_opened'](self)

    def _default_config(self):
        import rope.base.default_config
        import inspect
        return inspect.getsource(rope.base.default_config)

    def _init_other_parts(self):
        # Forcing the creation of `self.pycore` to register observers
        self.pycore

    def is_ignored(self, resource):
        return self.ignored.does_match(resource)

    def sync(self):
        """Closes project open resources"""
        self.close()

    def close(self):
        """Closes project open resources"""
        self.data_files.write()

    def set(self, key, value):
        """Set the `key` preference to `value`"""
        self.prefs.set(key, value)

    @property
    def ropefolder(self):
        if self._ropefolder_name is not None:
            return self.get_folder(self._ropefolder_name)

    def validate(self, folder=None):
        if folder is None:
            folder = self.root
        super(Project, self).validate(folder)

    root = property(lambda self: self.get_resource(''))
    address = property(lambda self: self._address)


class NoProject(_Project):
    """A null object for holding out of project files.

    This class is singleton use `get_no_project` global function
    """

    def __init__(self):
        fscommands = rope.base.fscommands.FileSystemCommands()
        super(NoProject, self).__init__(fscommands)

    def _get_resource_path(self, name):
        real_name = name.replace('/', os.path.sep)
        return _realpath(real_name)

    def get_resource(self, name):
        universal_name = _realpath(name).replace(os.path.sep, '/')
        return super(NoProject, self).get_resource(universal_name)

    def get_files(self):
        return []

    def get_python_files(self):
        return []

    _no_project = None


def get_no_project():
    if NoProject._no_project is None:
        NoProject._no_project = NoProject()
    return NoProject._no_project


class _FileListCacher(object):

    def __init__(self, project):
        self.project = project
        self.files = None
        rawobserver = resourceobserver.ResourceObserver(
            self._changed, self._invalid, self._invalid,
            self._invalid, self._invalid)
        self.project.add_observer(rawobserver)

    def get_files(self):
        if self.files is None:
            self.files = set()
            self._add_files(self.project.root)
        return self.files

    def _add_files(self, folder):
        for child in folder.get_children():
            if child.is_folder():
                self._add_files(child)
            elif not self.project.is_ignored(child):
                self.files.add(child)

    def _changed(self, resource):
        if resource.is_folder():
            self.files = None

    def _invalid(self, resource, new_resource=None):
        self.files = None


class _DataFiles(object):

    def __init__(self, project):
        self.project = project
        self.hooks = []

    def read_data(self, name, compress=False, import_=False):
        if self.project.ropefolder is None:
            return None
        compress = compress and self._can_compress()
        opener = self._get_opener(compress)
        file = self._get_file(name, compress)
        if not compress and import_:
            self._import_old_files(name)
        if file.exists():
            input = opener(file.real_path, 'rb')
            try:
                result = []
                try:
                    while True:
                        result.append(pickle.load(input))
                except EOFError:
                    pass
                if len(result) == 1:
                    return result[0]
                if len(result) > 1:
                    return result
            finally:
                input.close()

    def write_data(self, name, data, compress=False):
        if self.project.ropefolder is not None:
            compress = compress and self._can_compress()
            file = self._get_file(name, compress)
            opener = self._get_opener(compress)
            output = opener(file.real_path, 'wb')
            try:
                pickle.dump(data, output, 2)
            finally:
                output.close()

    def add_write_hook(self, hook):
        self.hooks.append(hook)

    def write(self):
        for hook in self.hooks:
            hook()

    def _can_compress(self):
        try:
            import gzip  # noqa
            return True
        except ImportError:
            return False

    def _import_old_files(self, name):
        old = self._get_file(name + '.pickle', False)
        new = self._get_file(name, False)
        if old.exists() and not new.exists():
            shutil.move(old.real_path, new.real_path)

    def _get_opener(self, compress):
        if compress:
            try:
                import gzip
                return gzip.open
            except ImportError:
                pass
        return open

    def _get_file(self, name, compress):
        path = self.project.ropefolder.path + '/' + name
        if compress:
            path += '.gz'
        return self.project.get_file(path)


def _realpath(path):
    """Return the real path of `path`

    Is equivalent to ``realpath(abspath(expanduser(path)))``.

    Of the particular notice is the hack dealing with the unfortunate
    sitaution of running native-Windows python (os.name == 'nt') inside
    of Cygwin (abspath starts with '/'), which apparently normal
    os.path.realpath completely messes up.

    """
    # there is a bug in cygwin for os.path.abspath() for abs paths
    if sys.platform == 'cygwin':
        if path[1:3] == ':\\':
            return path
        elif path[1:3] == ':/':
            path = "/cygdrive/" + path[0] + path[2:]
        return os.path.abspath(os.path.expanduser(path))
    return os.path.realpath(os.path.abspath(os.path.expanduser(path)))


def _find_module_in_folder(folder, modname):
    module = folder
    packages = modname.split('.')
    for pkg in packages[:-1]:
        if module.is_folder() and module.has_child(pkg):
            module = module.get_child(pkg)
        else:
            return None
    if module.is_folder():
        if module.has_child(packages[-1]) and \
           module.get_child(packages[-1]).is_folder():
            return module.get_child(packages[-1])
        elif module.has_child(packages[-1] + '.py') and \
                not module.get_child(packages[-1] + '.py').is_folder():
            return module.get_child(packages[-1] + '.py')
