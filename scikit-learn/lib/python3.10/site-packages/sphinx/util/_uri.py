from __future__ import annotations

from urllib.parse import parse_qsl, quote_plus, urlencode, urlsplit, urlunsplit


def encode_uri(uri: str) -> str:
    scheme, netloc, url, query, fragment = urlsplit(uri)
    netloc = netloc.encode('idna').decode('ascii')
    url = quote_plus(url.encode(), '/')
    query = urlencode([(q, v.encode()) for (q, v) in parse_qsl(query)])
    return urlunsplit((scheme, netloc, url, query, fragment))


def is_url(url: str) -> bool:
    """Check *url* is URL or not."""
    return bool(url) and '://' in url
