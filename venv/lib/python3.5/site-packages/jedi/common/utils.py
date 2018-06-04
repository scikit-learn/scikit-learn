import os


def traverse_parents(path, include_current=False):
    if not include_current:
        path = os.path.dirname(path)

    previous = None
    while previous != path:
        yield path
        previous = path
        path = os.path.dirname(path)
