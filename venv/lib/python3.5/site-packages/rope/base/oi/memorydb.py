from rope.base.oi import objectdb


class MemoryDB(objectdb.FileDict):

    def __init__(self, project, persist=None):
        self.project = project
        self._persist = persist
        self.files = self
        self._load_files()
        self.project.data_files.add_write_hook(self.write)

    def _load_files(self):
        self._files = {}
        if self.persist:
            result = self.project.data_files.read_data(
                'objectdb', compress=self.compress, import_=True)
            if result is not None:
                self._files = result

    def keys(self):
        return self._files.keys()

    def __iter__(self):
        for f in self._files:
            yield f

    def __len__(self):
        return len(self._files)

    def __setitem__(self):
        raise NotImplementedError()

    def __contains__(self, key):
        return key in self._files

    def __getitem__(self, key):
        return FileInfo(self._files[key])

    def create(self, path):
        self._files[path] = {}

    def rename(self, file, newfile):
        if file not in self._files:
            return
        self._files[newfile] = self._files[file]
        del self[file]

    def __delitem__(self, file):
        del self._files[file]

    def write(self):
        if self.persist:
            self.project.data_files.write_data('objectdb', self._files,
                                               self.compress)

    @property
    def compress(self):
        return self.project.prefs.get('compress_objectdb', False)

    @property
    def persist(self):
        if self._persist is not None:
            return self._persist
        else:
            return self.project.prefs.get('save_objectdb', False)


class FileInfo(objectdb.FileInfo):

    def __init__(self, scopes):
        self.scopes = scopes

    def create_scope(self, key):
        self.scopes[key] = ScopeInfo()

    def keys(self):
        return self.scopes.keys()

    def __contains__(self, key):
        return key in self.scopes

    def __getitem__(self, key):
        return self.scopes[key]

    def __delitem__(self, key):
        del self.scopes[key]

    def __iter__(self):
        for s in self.scopes:
            yield s

    def __len__(self):
        return len(self.scopes)

    def __setitem__(self):
        raise NotImplementedError()



class ScopeInfo(objectdb.ScopeInfo):

    def __init__(self):
        self.call_info = {}
        self.per_name = {}

    def get_per_name(self, name):
        return self.per_name.get(name, None)

    def save_per_name(self, name, value):
        self.per_name[name] = value

    def get_returned(self, parameters):
        return self.call_info.get(parameters, None)

    def get_call_infos(self):
        for args, returned in self.call_info.items():
            yield objectdb.CallInfo(args, returned)

    def add_call(self, parameters, returned):
        self.call_info[parameters] = returned

    def __getstate__(self):
        return (self.call_info, self.per_name)

    def __setstate__(self, data):
        self.call_info, self.per_name = data
