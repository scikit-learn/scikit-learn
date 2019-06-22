

anchor_attrib = '_yaml_anchor'


class Anchor(object):
    __slots__ = 'value', 'always_dump'
    attrib = anchor_attrib

    def __init__(self):
        # type: () -> None
        self.value = None
        self.always_dump = False
