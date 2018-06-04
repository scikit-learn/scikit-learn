from __future__ import print_function
try:
    from collections import MutableMapping
except ImportError:
    from UserDict import DictMixin as MutableMapping


class ObjectDB(object):

    def __init__(self, db, validation):
        self.db = db
        self.validation = validation
        self.observers = []
        self.files = db.files

    def validate_files(self):
        for file in list(self.files):
            if not self.validation.is_file_valid(file):
                del self.files[file]
                self._file_removed(file)

    def validate_file(self, file):
        if file not in self.files:
            return
        for key in list(self.files[file]):
            if not self.validation.is_scope_valid(file, key):
                del self.files[file][key]

    def file_moved(self, file, newfile):
        if file not in self.files:
            return
        self.files.rename(file, newfile)
        self._file_removed(file)
        self._file_added(newfile)

    def get_files(self):
        return self.files.keys()

    def get_returned(self, path, key, args):
        scope_info = self._get_scope_info(path, key, readonly=True)
        result = scope_info.get_returned(args)
        if self.validation.is_value_valid(result):
            return result

    def get_pername(self, path, key, name):
        scope_info = self._get_scope_info(path, key, readonly=True)
        result = scope_info.get_per_name(name)
        if self.validation.is_value_valid(result):
            return result

    def get_callinfos(self, path, key):
        scope_info = self._get_scope_info(path, key, readonly=True)
        return scope_info.get_call_infos()

    def add_callinfo(self, path, key, args, returned):
        scope_info = self._get_scope_info(path, key, readonly=False)
        old_returned = scope_info.get_returned(args)
        if self.validation.is_more_valid(returned, old_returned):
            scope_info.add_call(args, returned)

    def add_pername(self, path, key, name, value):
        scope_info = self._get_scope_info(path, key, readonly=False)
        old_value = scope_info.get_per_name(name)
        if self.validation.is_more_valid(value, old_value):
            scope_info.save_per_name(name, value)

    def add_file_list_observer(self, observer):
        self.observers.append(observer)

    def write(self):
        self.db.write()

    def _get_scope_info(self, path, key, readonly=True):
        if path not in self.files:
            if readonly:
                return _NullScopeInfo()
            self.files.create(path)
            self._file_added(path)
        if key not in self.files[path]:
            if readonly:
                return _NullScopeInfo()
            self.files[path].create_scope(key)
        result = self.files[path][key]
        if isinstance(result, dict):
            print(self.files, self.files[path], self.files[path][key])
        return result

    def _file_removed(self, path):
        for observer in self.observers:
            observer.removed(path)

    def _file_added(self, path):
        for observer in self.observers:
            observer.added(path)

    def __str__(self):
        scope_count = 0
        for file_dict in self.files.values():
            scope_count += len(file_dict)
        return 'ObjectDB holds %s file and %s scope infos' % \
               (len(self.files), scope_count)


class _NullScopeInfo(object):

    def __init__(self, error_on_write=True):
        self.error_on_write = error_on_write

    def get_per_name(self, name):
        pass

    def save_per_name(self, name, value):
        if self.error_on_write:
            raise NotImplementedError()

    def get_returned(self, parameters):
        pass

    def get_call_infos(self):
        return []

    def add_call(self, parameters, returned):
        if self.error_on_write:
            raise NotImplementedError()


class FileInfo(MutableMapping):

    def create_scope(self, key):
        pass


class FileDict(MutableMapping):

    def create(self, key):
        pass

    def rename(self, key, new_key):
        pass


class ScopeInfo(object):

    def get_per_name(self, name):
        pass

    def save_per_name(self, name, value):
        pass

    def get_returned(self, parameters):
        pass

    def get_call_infos(self):
        pass

    def add_call(self, parameters, returned):
        pass


class CallInfo(object):

    def __init__(self, args, returned):
        self.args = args
        self.returned = returned

    def get_parameters(self):
        return self.args

    def get_returned(self):
        return self.returned


class FileListObserver(object):

    def added(self, path):
        pass

    def removed(self, path):
        pass
