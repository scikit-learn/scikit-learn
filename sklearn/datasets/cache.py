# vendored from https://gitlab.com/KOLANICH/Cache.py
# Unlicense

import sqlite3
try:
    import ujson as json
except:
    import json

from collections import OrderedDict
import itertools
from types import FunctionType
import typing

class Transformer:
    __slots__ = ("id", "forth", "back")
    
    def __init__(self, id:str, forth:FunctionType, back:FunctionType):
        self.id = id
        self.forth = forth
        self.back = back
    
    def __repr__(self):
        return self.__class__.__name__ + "<" + self.id + ">"


class Compressor(Transformer):
    pass


class compressorsMeta(type):
    def __new__(cls, className:str, parents, attrs, *args, **kwargs):
        compressorsTemp = []
        
        try:
            import zopfli.zlib as zlib
        except ImportError:
            import zlib
        
        zlibCompressor = zlib.compressobj(level=zlib.Z_BEST_COMPRESSION, wbits=-15)
        zlibDecompressor = zlib.decompressobj(wbits=-15)
        compressorsTemp.append(Compressor("deflate", zlibCompressor.compress, zlibDecompressor.decompress))
        
        try:
            import lz4.frame
            
            def lz4Compress(data:bytes):
                return lz4.frame.compress(
                    data, compression_level=lz4.frame.COMPRESSIONLEVEL_MAX,
                    block_size=lz4.frame.BLOCKSIZE_MAX4MB,
                    return_bytearray=True
                )
            
            def lz4Decompress(data:bytes):
                return lz4.frame.decompress(data, return_bytearray=True)
            
            compressorsTemp.append(Compressor("lz4", lz4Compress, lz4Decompress))
        except:
            pass
        
        try:
            import brotli
            compressorsTemp.append(Compressor("brotli", brotli.compress, brotli.decompress))
        except:
            pass
        
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
            
            def lzmaCompressFunc(data:bytes):
                return lzma.compress(data, **lzmaArgs)
            
            def lzmaDecompressFunc(data:bytes):
                return lzma.decompress(data, **lzmaArgs)
            
            compressorsTemp.append(Compressor("lzma", lzmaCompressFunc, lzmaDecompressFunc))
        except:
            pass
        
        attrsNew = OrderedDict()
        attrsNew["none"] = None
        
        attrsNew.update((c.id, c) for c in compressorsTemp)
        
        attrsNew["_BEST"] = compressorsTemp[-1]
        attrsNew.update(attrs)
        
        return super().__new__(cls, className, parents, attrsNew, *args, **kwargs)


class compressors(metaclass=compressorsMeta):
    pass


class CacheMeta(type):
    def __new__(cls, className:str, parents, attrs, *args, **kwargs):
        assert(len(parents) <= 1)
        if "_appendSerializers" in attrs:
            if parents:
                attrs["serializers"] = list(parents[0].serializers)
            else:
                attrs["serializers"] = tuple()
            
            attrs["serializers"].extend(attrs["_appendSerializers"])
            attrs["_serializersChainSeqId"] = tuple((s.id for s in attrs["serializers"]))
        
        return super().__new__(cls, className, parents, attrs, *args, **kwargs)


class BlobCache(metaclass=CacheMeta):
    """Just a simple SQLite-based cache"""
    __slots__ = ("compressor", "db", "path")
    TABLE_NAME = "cache"
    META_TABLE_NAME = "metadata"
    serializers = ()
    _serializersChainSeqId = ()
    
    def __init__(self, path:typing.Union["pathlib.Path", str]="./cache", compressor:typing.Optional[typing.Union[str, Compressor]]=None) -> None:
        if isinstance(compressor, str):
            compressor = getattr(compressors, compressor)
        elif compressor is True:
            compressor = compressors._BEST
        
        self.path = path
        self.db = None
        self.compressor = compressor
    
    def createMetaDataTable(self) -> None:
        self.db.executescript(
            r"""create table `""" + self.__class__.META_TABLE_NAME + """` (
                key TEXT PRIMARY KEY,
                val TEXT
            );
            """
        )
    
    def createDataTable(self) -> None:
        self.db.executescript(
            r"""create table `""" + self.__class__.TABLE_NAME + r"""` (
                key TEXT PRIMARY KEY,
                val BLOB
            );
            """
        )
    
    def __len__(self) -> int:
        cur = self.db.execute("select count(*) from `" + self.__class__.TABLE_NAME + "`;")
        res = next(cur)[0]
        cur.close()
        return res
    
    def __contains__(self, key:str) -> bool:
        return self[key] is not None
    
    def getMetadata(self, key:str) -> str:
        return next(
            self.db.execute(
                "select `val` from `" + self.__class__.META_TABLE_NAME + "` where `key` = ?;",
                (key,)
            )
        )[0]
    
    def setMetadata(self, key:str, value:str) -> None:
        self.db.execute(
            "insert or replace into `" + self.__class__.META_TABLE_NAME + "` (`key`, `val`) values (?, ?);",
            (key, value)
        )
        self.db.commit()
    
    def __getitem__(self, key:str) -> object:
        cur = self.db.execute("select `val` from `" + self.__class__.TABLE_NAME + "` where `key` = ?;", (key,))
        res = tuple(cur)
        if not res:
            return None
        
        res = res[0][0]
        
        transformers = self.__class__.serializers
        if self.compressor:
            transformers = itertools.chain((self.compressor,), transformers)
        
        for serializer in transformers:
            prevRes = res
            #sys.stdout.write("?"+"<= "+repr(res)+"<= "+serializer.id)
            res = serializer.back(res)
            # sys.stdout.write("\b")
            #print(repr(res), "<=", repr(prevRes))
        # print("\n")
        
        cur.close()
        return res
    
    def __setitem__(self, key:str, val:object) -> None:
        if val is None:
            del(self[key])
        else:
            transformers = reversed(self.__class__.serializers)
            if self.compressor:
                transformers = itertools.chain(
                    transformers, (self.compressor,))
            
            for serializer in transformers:
                #sys.stdout.write(repr(val)+" => "+serializer.id)
                val = serializer.forth(val)
                #sys.stdout.write(" => "+repr(val)+"\n")
            #print("\n")
            
            return self.db.execute(
                "insert or replace into `" + self.__class__.TABLE_NAME + "` (`key`, `val`) values (?, ?);",
                (key, val)
            ) and self.db.commit()
    
    def __delitem__(self, key) -> None:
        self.db.execute(
            "delete from `" + self.__class__.TABLE_NAME + "` where `key` = ?;",
            (key,)
        ) and self.db.commit()
    
    def isInitialized(self) -> bool:
        return next(self.db.execute("SELECT count(*) FROM `sqlite_master` WHERE `type`='table' AND `name`='" + self.__class__.META_TABLE_NAME + "';"))[0]
    
    def __enter__(self) -> "BlobCache":
        self.db = sqlite3.connect(str(self.path))
        
        if not self.isInitialized():
            self.createMetaDataTable()
            self.setMetadata("compression", self.compressor.id)
            self.setMetadata("serializers", json.dumps(self.__class__._serializersChainSeqId))
            self.createDataTable()
        else:
            compressionId = self.getMetadata("compression")
            serializersTemp = tuple(json.loads(self.getMetadata("serializers")))
            if serializersTemp != self.__class__._serializersChainSeqId:
                raise ValueError("This DB serializers chain doesn't match the chain used in the class: " + repr(serializersTemp),)
            
            self.compressor = getattr(compressors, compressionId)
        
        return self
    
    def __exit__(self, exc_class, exc, traceback) -> None:
        self.db.commit()
        self.db.close()
        self.db = None
    
    def __del__(self) -> None:
        if self.db is not None:
            self.__exit__(None, None, None)
    
    def empty(self) -> None:
        """Empties the DB"""
        self.db.execute("drop table `" + self.__class__.TABLE_NAME + "`;")
        self.createDataTable()
        self.db.commit()
    
    def vacuum(self) -> None:
        self.db.execute("reindex;")
        self.db.execute("vacuum;")


class StringCache(BlobCache):
    _appendSerializers = (Transformer("utf-8", lambda d: d.encode("utf-8"), lambda d: d.decode("utf-8")),)


class JSONCache(StringCache):
    _appendSerializers = (Transformer("json", lambda d: json.dumps(d), lambda d: json.loads(d)),)


Cache = JSONCache

try:
    import bson.json_util
    
    class BSONCache(BlobCache):
        _appendSerializers = (Transformer("bson", lambda d: bson.json_util.dumps(d), lambda d: bson.json_util.loads(d)),)
    Cache = BSONCache
except ImportError:
    pass

try:
    import msgpack
    msgpackPacker = msgpack.Packer(use_bin_type=True, strict_types=True)
    
    class MSGPackCache(BlobCache):
        _appendSerializers = (Transformer("msgpack", lambda d: msgpackPacker.pack(d), lambda d: msgpack.unpackb(d, raw=False)),)
    Cache = MSGPackCache
except ImportError:
    pass
