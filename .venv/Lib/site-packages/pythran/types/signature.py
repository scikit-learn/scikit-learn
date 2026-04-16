from pythran.typing import List, Dict, Set, Fun, TypeVar, Type
from pythran.typing import Union, Iterable


def type_dependencies(t):

    if isinstance(t, TypeVar):
        return {t}
    else:
        return set().union(*[type_dependencies(arg)
                             for arg in getattr(t, '__args__', [])])


def dep_builder(type_var, ppal_index, index, t, self, node):
    if isinstance(t, TypeVar):
        if t is type_var:
            # FIXME: this is the second part of the hack below,
            # FIXME: there's no reason why self.result[node.args[index]]
            # FIXME: would still be valid in translated context
            return lambda arg: (arg
                                if index == ppal_index
                                else self.result[node.args[index]])
    elif isinstance(t, (List, Set, Iterable, Dict, Type)):
        return lambda arg: self.builder.IteratorContentType(
            dep_builder(type_var,
                        ppal_index,
                        index,
                        t.__args__[0],
                        self, node)(arg))
    assert False, t


class InfeasibleCombiner(Exception):
    pass


def path_to(self, t, deps_builders, node):

    if isinstance(t, TypeVar):
        if t in deps_builders:
            return deps_builders[t]
        else:
            raise InfeasibleCombiner()
    if isinstance(t, Type):
        return lambda arg: self.builder.TypeType(
            path_to(self, t.__args__[0], deps_builders, node)(arg))
    if isinstance(t, List):
        return lambda arg: self.builder.ListType(
            path_to(self, t.__args__[0], deps_builders, node)(arg))
    if isinstance(t, Set):
        return lambda arg: self.builder.SetType(
            path_to(self, t.__args__[0], deps_builders, node)(arg))
    if isinstance(t, Dict):
        return lambda arg: self.builder.DictType(
            path_to(self, t.__args__[0], deps_builders, node)(arg),
            path_to(self, t.__args__[1], deps_builders, node)(arg),
        )
    if isinstance(t, Fun):
        raise InfeasibleCombiner()
    if isinstance(t, Iterable):  # FIXME?
        raise InfeasibleCombiner()

    assert False, (t, t.mro())


def build_unary_op(deps, args, self, node):
    # FIXME: this is a hack, because only the fist dep gets translated
    # FIXME: in case of interprocedural translation
    # FIXME: this was the original behavior...
    ppal_index = sorted(deps.values())[0][0][0]
    deps_builders = {dep: dep_builder(dep,
                                      ppal_index,
                                      *src[0],
                                      self=self,
                                      node=node)
                     for dep, src in deps.items()}
    return path_to(self, args[0], deps_builders, node), ppal_index


def build_combiner(signature, deps):
    sig_args = signature.__args__[:-1]

    def combiner(self, node):
        if deps and len(node.args) == len(sig_args):
            try:

                unary_op, main_index = build_unary_op(deps,
                                                      sig_args,
                                                      self, node)

                self.combine(
                    node.args[0],
                    unary_op,
                    node.args[main_index])
            except InfeasibleCombiner:
                pass

    return combiner


def extract_combiner(signature):
    if not isinstance(signature, (Fun, Union)):
        return None

    if type(signature) is Union:
        combiners = [extract_combiner(up)
                     for up in signature.__args__]
        combiners = [cb for cb in combiners if cb]

        def combiner(self, node):
            for cb in combiners:
                cb(self, node)
        return combiner

    args = signature.__args__[:-1]

    if not args:
        return None

    deps = type_dependencies(args[0])

    if not deps:
        return None

    deps_src = dict()

    for i, arg in enumerate(args[1:]):
        arg_deps = type_dependencies(arg)
        common_deps = deps.intersection(arg_deps)
        for common_dep in common_deps:
            deps_src.setdefault(common_dep, []).append((i + 1, arg))

    return build_combiner(signature, deps_src)
