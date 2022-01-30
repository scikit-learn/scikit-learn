import sys

class _Method: ...

METHOD_CRYPT: _Method
METHOD_MD5: _Method
METHOD_SHA256: _Method
METHOD_SHA512: _Method
if sys.version_info >= (3, 7):
    METHOD_BLOWFISH: _Method

methods: list[_Method]

if sys.version_info >= (3, 7):
    def mksalt(method: _Method | None = ..., *, rounds: int | None = ...) -> str: ...

else:
    def mksalt(method: _Method | None = ...) -> str: ...

def crypt(word: str, salt: str | _Method | None = ...) -> str: ...
