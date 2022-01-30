from typing import Optional


class StubInfo:
    def __init__(self, name: str, py_version: Optional[int] = None) -> None:
        self.name = name
        # If None, compatible with py2+py3, if 2/3, only compatible with py2/py3
        self.py_version = py_version


def is_legacy_bundled_package(prefix: str, py_version: int) -> bool:
    if prefix not in legacy_bundled_packages:
        return False
    package_ver = legacy_bundled_packages[prefix].py_version
    return package_ver is None or package_ver == py_version


# Stubs for these third-party packages used to be shipped with mypy.
#
# Map package name to PyPI stub distribution name.
#
# Package name can have one or two components ('a' or 'a.b').
legacy_bundled_packages = {
    'aiofiles': StubInfo('types-aiofiles', py_version=3),
    'atomicwrites': StubInfo('types-atomicwrites'),
    'attr': StubInfo('types-attrs'),
    'backports': StubInfo('types-backports'),
    'backports_abc': StubInfo('types-backports_abc'),
    'bleach': StubInfo('types-bleach'),
    'boto': StubInfo('types-boto'),
    'cachetools': StubInfo('types-cachetools'),
    'certifi': StubInfo('types-certifi'),
    'characteristic': StubInfo('types-characteristic'),
    'chardet': StubInfo('types-chardet'),
    'click_spinner': StubInfo('types-click-spinner'),
    'concurrent': StubInfo('types-futures', py_version=2),
    'contextvars': StubInfo('types-contextvars', py_version=3),
    'croniter': StubInfo('types-croniter'),
    'dataclasses': StubInfo('types-dataclasses', py_version=3),
    'dateparser': StubInfo('types-dateparser'),
    'datetimerange': StubInfo('types-DateTimeRange'),
    'dateutil': StubInfo('types-python-dateutil'),
    'decorator': StubInfo('types-decorator'),
    'deprecated': StubInfo('types-Deprecated'),
    'docutils': StubInfo('types-docutils', py_version=3),
    'emoji': StubInfo('types-emoji'),
    'enum': StubInfo('types-enum34', py_version=2),
    'fb303': StubInfo('types-fb303', py_version=2),
    'filelock': StubInfo('types-filelock', py_version=3),
    'first': StubInfo('types-first'),
    'freezegun': StubInfo('types-freezegun', py_version=3),
    'frozendict': StubInfo('types-frozendict', py_version=3),
    'geoip2': StubInfo('types-geoip2'),
    'gflags': StubInfo('types-python-gflags'),
    'google.protobuf': StubInfo('types-protobuf'),
    'ipaddress': StubInfo('types-ipaddress', py_version=2),
    'kazoo': StubInfo('types-kazoo', py_version=2),
    'markdown': StubInfo('types-Markdown'),
    'maxminddb': StubInfo('types-maxminddb'),
    'mock': StubInfo('types-mock'),
    'OpenSSL': StubInfo('types-pyOpenSSL'),
    'orjson': StubInfo('types-orjson', py_version=3),
    'paramiko': StubInfo('types-paramiko'),
    'pathlib2': StubInfo('types-pathlib2', py_version=2),
    'pkg_resources': StubInfo('types-setuptools', py_version=3),
    'polib': StubInfo('types-polib'),
    'pycurl': StubInfo('types-pycurl'),
    'pymssql': StubInfo('types-pymssql', py_version=2),
    'pymysql': StubInfo('types-PyMySQL'),
    'pyrfc3339': StubInfo('types-pyRFC3339', py_version=3),
    'python2': StubInfo('types-six'),
    'pytz': StubInfo('types-pytz'),
    'pyVmomi': StubInfo('types-pyvmomi'),
    'redis': StubInfo('types-redis'),
    'requests': StubInfo('types-requests'),
    'retry': StubInfo('types-retry'),
    'routes': StubInfo('types-Routes', py_version=2),
    'scribe': StubInfo('types-scribe', py_version=2),
    'simplejson': StubInfo('types-simplejson'),
    'singledispatch': StubInfo('types-singledispatch'),
    'six': StubInfo('types-six'),
    'slugify': StubInfo('types-python-slugify'),
    'tabulate': StubInfo('types-tabulate'),
    'termcolor': StubInfo('types-termcolor'),
    'toml': StubInfo('types-toml'),
    'tornado': StubInfo('types-tornado', py_version=2),
    'typed_ast': StubInfo('types-typed-ast', py_version=3),
    'tzlocal': StubInfo('types-tzlocal'),
    'ujson': StubInfo('types-ujson'),
    'waitress': StubInfo('types-waitress', py_version=3),
    'yaml': StubInfo('types-PyYAML'),
}
