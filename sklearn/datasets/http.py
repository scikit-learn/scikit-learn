from . import cache

headers = {
    #"User-Agent": "scikit-learn datasets fetcher",
    "Connection": "keep-alive",
}

try:
    import requests

    def setupSession():
        reqSess = requests.Session()

        try:
            import hyper.http20.response
            headers["Accept-Encoding"] = [k.decode("utf-8") for k in hyper.http20.response.decompressors]

            if b"br" not in hyper.http20.response.decompressors:
                try:
                    import brotli
                    hyper.http20.response.decompressors[b"br"] = brotli.Decompressor
                except:
                    pass
            if b"br" in hyper.http20.response.decompressors:
                headers["Accept-Encoding"].append("br")

            from hyper.contrib import HTTP20Adapter
            ad = HTTP20Adapter()
            reqSess.mount("https://", ad)
            reqSess.mount("http://", ad)
            hyper.h2.settings.ENABLE_PUSH = False
        except:
            pass

        headers["Accept-Encoding"] = ", ".join(headers["Accept-Encoding"])
        reqSess.headers.update(headers)
        return reqSess

    reqSess = setupSession()

    def doGet(uri, binary=False):
        res = reqSess.get(uri)
        if binary:
            return res.raw
        return res.text
except:
    try:
        # Python 3+
        from urllib.request import urlopen, Request
    except ImportError:
        # Python 2
        from urllib2 import urlopen, Request

    headers["Accept-Encoding"] = ["gzip", "deflate"]
    import gzip
    import zlib

    decompressors = {
        "gzip": gzip.decompress,
        "deflate": zlib.decompress
    }

    try:
        import brotli
        decompressors["br"] = brotli.Decompressor
    except:
        pass

    headers["Accept-Encoding"] = ", ".join(headers["Accept-Encoding"])

    def doGet(uri, binary=False):
        req = Request(uri)

        for k, v in headers.items():
            req.add_header(k, v)

        res = urlopen(req)
        enc = res.getheader('Content-Encoding')
        res = res.read()
        if enc in decompressors:
            res = decompressors[enc](res)

        if not binary:
            ct = a.getheader('Content-Type').split(";")
            encoding = "utf-8"
            for p in ct:
                splitted = p.split("=")
                if len(splitted) == 2:
                    k, v = splitted
                    k = k.strip()
                    v = v.strip()
                    if k == "charset":
                        encoding = v
            res = res.decode(encoding)
        return res


def get(uri, *args, force=False, binary=False, cacheFile="./scikit-learn_datasets_HTTP_cache.sqlite", **kwargs):
    with cache.StringCache(cacheFile, True) as hcache:
        res = hcache[uri]
        if res and not force:
            return res
        else:
            res = doGet(uri, binary)
            if res:
                hcache[uri] = res
            return res
