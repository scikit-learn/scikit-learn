"""Interfaces for accessing metadata.

We provide two implementations.
 * The "classic" file system implementation, which uses a directory
   structure of files.
 * A hokey sqlite backed implementation, which basically simulates
   the file system in an effort to work around poor file system performance
   on OS X.
"""

import binascii
import os
import time

from abc import abstractmethod
from typing import List, Iterable, Any, Optional
from typing_extensions import TYPE_CHECKING
if TYPE_CHECKING:
    # We avoid importing sqlite3 unless we are using it so we can mostly work
    # on semi-broken pythons that are missing it.
    import sqlite3


class MetadataStore:
    """Generic interface for metadata storage."""

    @abstractmethod
    def getmtime(self, name: str) -> float:
        """Read the mtime of a metadata entry..

        Raises FileNotFound if the entry does not exist.
        """
        pass

    @abstractmethod
    def read(self, name: str) -> str:
        """Read the contents of a metadata entry.

        Raises FileNotFound if the entry does not exist.
        """
        pass

    @abstractmethod
    def write(self, name: str, data: str, mtime: Optional[float] = None) -> bool:
        """Write a metadata entry.

        If mtime is specified, set it as the mtime of the entry. Otherwise,
        the current time is used.

        Returns True if the entry is successfully written, False otherwise.
        """

    @abstractmethod
    def remove(self, name: str) -> None:
        """Delete a metadata entry"""
        pass

    @abstractmethod
    def commit(self) -> None:
        """If the backing store requires a commit, do it.

        But N.B. that this is not *guaranteed* to do anything, and
        there is no guarantee that changes are not made until it is
        called.
        """
        pass

    @abstractmethod
    def list_all(self) -> Iterable[str]: ...


def random_string() -> str:
    return binascii.hexlify(os.urandom(8)).decode('ascii')


class FilesystemMetadataStore(MetadataStore):
    def __init__(self, cache_dir_prefix: str) -> None:
        # We check startswith instead of equality because the version
        # will have already been appended by the time the cache dir is
        # passed here.
        if cache_dir_prefix.startswith(os.devnull):
            self.cache_dir_prefix = None
        else:
            self.cache_dir_prefix = cache_dir_prefix

    def getmtime(self, name: str) -> float:
        if not self.cache_dir_prefix:
            raise FileNotFoundError()

        return int(os.path.getmtime(os.path.join(self.cache_dir_prefix, name)))

    def read(self, name: str) -> str:
        assert os.path.normpath(name) != os.path.abspath(name), "Don't use absolute paths!"

        if not self.cache_dir_prefix:
            raise FileNotFoundError()

        with open(os.path.join(self.cache_dir_prefix, name), 'r') as f:
            return f.read()

    def write(self, name: str, data: str, mtime: Optional[float] = None) -> bool:
        assert os.path.normpath(name) != os.path.abspath(name), "Don't use absolute paths!"

        if not self.cache_dir_prefix:
            return False

        path = os.path.join(self.cache_dir_prefix, name)
        tmp_filename = path + '.' + random_string()
        try:
            os.makedirs(os.path.dirname(path), exist_ok=True)
            with open(tmp_filename, 'w') as f:
                f.write(data)
            os.replace(tmp_filename, path)
            if mtime is not None:
                os.utime(path, times=(mtime, mtime))

        except os.error:
            return False
        return True

    def remove(self, name: str) -> None:
        if not self.cache_dir_prefix:
            raise FileNotFoundError()

        os.remove(os.path.join(self.cache_dir_prefix, name))

    def commit(self) -> None:
        pass

    def list_all(self) -> Iterable[str]:
        if not self.cache_dir_prefix:
            return

        for dir, _, files in os.walk(self.cache_dir_prefix):
            dir = os.path.relpath(dir, self.cache_dir_prefix)
            for file in files:
                yield os.path.join(dir, file)


SCHEMA = '''
CREATE TABLE IF NOT EXISTS files (
    path TEXT UNIQUE NOT NULL,
    mtime REAL,
    data TEXT
);
CREATE INDEX IF NOT EXISTS path_idx on files(path);
'''
# No migrations yet
MIGRATIONS: List[str] = []


def connect_db(db_file: str) -> 'sqlite3.Connection':
    import sqlite3.dbapi2

    db = sqlite3.dbapi2.connect(db_file)
    db.executescript(SCHEMA)
    for migr in MIGRATIONS:
        try:
            db.executescript(migr)
        except sqlite3.OperationalError:
            pass
    return db


class SqliteMetadataStore(MetadataStore):
    def __init__(self, cache_dir_prefix: str) -> None:
        # We check startswith instead of equality because the version
        # will have already been appended by the time the cache dir is
        # passed here.
        if cache_dir_prefix.startswith(os.devnull):
            self.db = None
            return

        os.makedirs(cache_dir_prefix, exist_ok=True)
        self.db = connect_db(os.path.join(cache_dir_prefix, 'cache.db'))

    def _query(self, name: str, field: str) -> Any:
        # Raises FileNotFound for consistency with the file system version
        if not self.db:
            raise FileNotFoundError()

        cur = self.db.execute('SELECT {} FROM files WHERE path = ?'.format(field), (name,))
        results = cur.fetchall()
        if not results:
            raise FileNotFoundError()
        assert len(results) == 1
        return results[0][0]

    def getmtime(self, name: str) -> float:
        return self._query(name, 'mtime')

    def read(self, name: str) -> str:
        return self._query(name, 'data')

    def write(self, name: str, data: str, mtime: Optional[float] = None) -> bool:
        import sqlite3

        if not self.db:
            return False
        try:
            if mtime is None:
                mtime = time.time()
            self.db.execute('INSERT OR REPLACE INTO files(path, mtime, data) VALUES(?, ?, ?)',
                            (name, mtime, data))
        except sqlite3.OperationalError:
            return False
        return True

    def remove(self, name: str) -> None:
        if not self.db:
            raise FileNotFoundError()

        self.db.execute('DELETE FROM files WHERE path = ?', (name,))

    def commit(self) -> None:
        if self.db:
            self.db.commit()

    def list_all(self) -> Iterable[str]:
        if self.db:
            for row in self.db.execute('SELECT path FROM files'):
                yield row[0]
