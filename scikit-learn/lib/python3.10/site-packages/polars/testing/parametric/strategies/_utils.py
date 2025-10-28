from typing import Any


def flexhash(elem: Any) -> int:
    """
    Hashing function that also handles lists and dictionaries.

    Used for `unique` check in nested strategies.
    """
    if isinstance(elem, list):
        return hash(tuple(flexhash(e) for e in elem))
    elif isinstance(elem, dict):
        return hash(tuple((k, flexhash(v)) for k, v in elem.items()))
    return hash(elem)
