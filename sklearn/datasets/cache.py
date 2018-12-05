# vendored from https://gitlab.com/KOLANICH/Cache.py
# Unlicense

import sqlite3
try:
    import ujson as json
except:
    import json

from collections import OrderedDict

class Transformer:
    __slots__ = ("id", "forth", "back")
    
    def __init__(self, id, forth, back):
        self.id = id
        self.forth = forth
        self.back = back
    
    def __repr__(self):
        return self.__class__.__name__+"<"+self.id+">"

class Compressor(Transformer):
    pass

class compressorsMeta(type):
    def __new__(cls, className, parents, attrs, *args, **kwargs):
        def dummy(x):
            return x
        compressorsTemp = [Compressor("none", dummy, dummy)]
        
        import gzip
        compressorsTemp.append(Compressor("gzip", gzip.compress, gzip.decompress))
        
        try:
            try:
                import lzma
            except ImportError:
                import backports.lzma
            
            lzmaArgs = {
                "format": lzma.FORMAT_RAW,
                "filters": (
                    {"id": lzma.FILTER_DELTA, "dist": 5},
                    {"id": lzma.FILTER_LZMA2, "preset": 9 | lzma.PRESET_EXTREME}
                )
            }
            
            def compressFunc(data):
                return lzma.compress(data, **lzmaArgs)
            
            def decompressFunc(data):
                return lzma.decompress(data, **lzmaArgs)
            
            compressorsTemp.append(Compressor("lzma", compressFunc, decompressFunc))
        except:
            pass
        
        attrsNew = OrderedDict(((c.id, c) for c in compressorsTemp))
        attrsNew["_BEST"] = compressorsTemp[-1]
        attrsNew.update(attrs)
        
        return super().__new__(cls, className, parents, attrsNew, *args, **kwargs)


class compressors(metaclass=compressorsMeta):
    pass

class CacheMeta(type):
    def __new__(cls, className, parents, attrs, *args, **kwargs):
        assert(len(parents) <= 1)
        if "_appendSerializers" in attrs:
            if parents:
                attrs["serializers"]=list(parents[0].serializers)
            else:
                attrs["serializers"]=tuple()
            
            attrs["serializers"].extend(attrs["_appendSerializers"])
            attrs["_serializersChainSeqId"]=tuple((s.id for s in attrs["serializers"]))
        
        return super().__new__(cls, className, parents, attrs, *args, **kwargs)

class BlobCache(metaclass=CacheMeta):
    """Just a simple SQLite-based cache"""
    __slots__ = ("compressor", "db", "path")
    TABLE_NAME = "cache"
    META_TABLE_NAME = "metadata"
    serializers=()
    
    def __init__(self, path="./cache", compressor=None):
        if isinstance(compressor, str):
            compressor = getattr(compressors, compressor)
        elif compressor is True:
            compressor = compressors._BEST
        elif compressor is None:
            compressor = compressors.none
        
        self.path=path
        self.db=None
        self.compressor = compressor

    def createMetaDataTable(self):
        self.db.executescript(
            r"""create table `"""+self.__class__.META_TABLE_NAME+"""` (
                key TEXT PRIMARY KEY,
                val TEXT
            );
            """
        )
    
    def createDataTable(self):
        self.db.executescript(
            r"""create table `""" + self.__class__.TABLE_NAME + r"""` (
                key TEXT PRIMARY KEY,
                val BLOB
            );
            """
        )

    def __len__(self):
        cur = self.db.execute("select count(*) from `" + self.__class__.TABLE_NAME + "`;")
        res = next(cur)[0]
        cur.close()
        return res

    def __contains__(self, key):
        return self[key] != None

    def getMetadata(self, key):
        return next(self.db.execute(
            "select `val` from `" + self.__class__.META_TABLE_NAME + "` where `key` = ?;",
            (key,)
        ))[0]

    def setMetadata(self, key, value):
        self.db.execute(
            "insert or replace into `"+self.__class__.META_TABLE_NAME + "` (`key`, `val`) values (?, ?);",
            (key, value)
        )
        self.db.commit()

    def __getitem__(self, key):
        cur = self.db.execute(
            "select `val` from `"+self.__class__.TABLE_NAME+"` where `key` = ?;",
            (key,)
        )
        res = tuple(cur)
        try:
            res = res[0][0]
            res = self.compressor.back(res)
            for serializer in reversed(self.__class__.serializers):
                res = serializer.back(res)
        except:
            res = None
        cur.close()
        return res

    def __setitem__(self, key, val):
        if val is None:
            del(self[key])
        else:
            for serializer in self.__class__.serializers:
                val = serializer.forth(val)
            
            return self.db.execute(
                "insert or replace into `" + self.__class__.TABLE_NAME + "` (`key`, `val`) values (?, ?);",
                (key, self.compressor.forth(val))
            ) and self.db.commit()

    def __delitem__(self, key):
        return self.db.execute(
            "delete from `" + self.__class__.TABLE_NAME + "` where `key` = ?;",
            (key,)
        ) and self.db.commit()
    
    def isInitialized(self):
        return next(self.db.execute("SELECT count(*) FROM `sqlite_master` WHERE `type`='table' AND `name`='"+self.__class__.META_TABLE_NAME+"';"))[0]
    
    def __enter__(self):
        self.db = sqlite3.connect(self.path)
        
        if not self.isInitialized():
            self.createMetaDataTable()
            self.setMetadata("compression", self.compressor.id)
            self.setMetadata("serializers", json.dumps(self.__class__._serializersChainSeqId))
            self.createDataTable()
        else:
            compressionId = self.getMetadata("compression")
            serializersTemp = tuple(json.loads(self.getMetadata("serializers")))
            if serializersTemp != self.__class__._serializersChainSeqId:
                raise ValueError("This DB serializers chain doesn't match the chain used in the class: " + repr(serializersTemp), )
            
            self.compressor = getattr(compressors, compressionId)
        
        return self

    def __exit__(self, exc_class, exc, traceback):
        self.db.commit()
        self.db.close()
        self.db = None
    
    def __del__(self):
        if self.db is not None:
            self.__exit__(None, None, None)

    def empty(self):
        """Empties the DB"""
        self.db.executescript("drop table `" + self.__class__.TABLE_NAME+"`;")
        self.createDataTable()
        self.db.commit()
        


class StringCache(BlobCache):
    _appendSerializers=(Transformer("utf-8", lambda d: d.encode("utf-8"), lambda d: d.decode("utf-8")),)

class JSONCache(StringCache):
    _appendSerializers=(Transformer("json", lambda d: json.dumps(d), lambda d: json.loads(d)),)
Cache = JSONCache

try:
    import bson
    class BSONCache(BlobCache):
        _appendSerializers=(Transformer("bson", lambda d: bson.dumps(d), lambda d: bson.loads(d)),)
    Cache = BSONCache
except ImportError:
    pass

try:
    import msgpack
    class MSGPackCache(BlobCache):
        _appendSerializers=(Transformer("msgpack", lambda d: msgpack.dumps(d), lambda d: msgpack.loads(d)),)
    Cache = MSGPackCache
except ImportError:
    pass
