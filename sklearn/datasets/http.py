from . import cache

headers = {
    #"User-Agent": "scikit-learn datasets fetcher",
    "Connection": "keep-alive",
}

import requests

def setup_session():
    requests_session = requests.Session()

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
        requests_session.mount("https://", ad)
        requests_session.mount("http://", ad)
        hyper.h2.settings.ENABLE_PUSH = False
    except:
        pass

    headers["Accept-Encoding"] = ", ".join(headers["Accept-Encoding"])
    requests_session.headers.update(headers)
    return requests_session

requests_session = setup_session()

def do_get(uri, binary=False):
    res = requests_session.get(uri)
    if binary:
        return res.raw
    return res.text


def get(uri, *args, force=False, binary=False, cacheFile="./scikit-learn_datasets_HTTP_cache.sqlite", **kwargs):
    with cache.StringCache(cacheFile, True) as hcache:
        res = hcache[uri]
        if res and not force:
            return res
        else:
            res = do_get(uri, binary)
            if res:
                hcache[uri] = res
            return res
