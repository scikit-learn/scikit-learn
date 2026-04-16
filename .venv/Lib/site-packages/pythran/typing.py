class FunMeta(type):

    def __getitem__(cls, item):
        return Fun(tuple(item[0]) + (item[1],))


class DictMeta(type):

    def __getitem__(cls, item):
        return Dict(item)


class UnionMeta(type):

    def __getitem__(cls, item):
        return Union(item)


class SetMeta(type):

    def __getitem__(cls, item):
        return Set(item)


class ListMeta(type):

    def __getitem__(cls, item):
        return List(item)


class TypeMeta(type):

    def __getitem__(cls, item):
        return Type(item)


class IterableMeta(type):

    def __getitem__(cls, item):
        return Iterable(item)


class GeneratorMeta(type):

    def __getitem__(cls, item):
        return Generator(item)


class TupleMeta(type):

    def __getitem__(cls, item):
        return Tuple(item)


class OptionalMeta(type):

    def __getitem__(cls, item):
        return Optional(item)


class NDArrayMeta(type):

    def __getitem__(cls, item):
        return NDArray(item)


class PointerMeta(type):

    def __getitem__(cls, item):
        return Pointer(item)


class TypeBase(type):

    def __new__(cls, args):
        return type.__new__(
            cls,
            cls.__name__,
            (object,),
            {'__args__': args if isinstance(args, tuple) else (args,)}
        )

    def __init__(self, *args, **kwargs):
        pass


class Fun(TypeBase, metaclass=FunMeta):
    pass


class Dict(TypeBase, metaclass=DictMeta):
    pass


class Union(TypeBase, metaclass=UnionMeta):
    pass


class Set(TypeBase, metaclass=SetMeta):
    pass


class List(TypeBase, metaclass=ListMeta):
    pass


class Iterable(TypeBase, metaclass=IterableMeta):
    pass


class Generator(TypeBase, metaclass=GeneratorMeta):
    pass


class Tuple(TypeBase, metaclass=TupleMeta):
    pass


class Optional(TypeBase, metaclass=OptionalMeta):
    pass


class NDArray(TypeBase, metaclass=NDArrayMeta):
    pass


class Pointer(TypeBase, metaclass=PointerMeta):
    pass


class TypeVar(object):

    def __init__(self, name):
        self.__name__ = name


class Sized(object):
    pass


class Any(object):
    pass


class File(object):
    pass


class Type(TypeBase, metaclass=TypeMeta):
    pass
