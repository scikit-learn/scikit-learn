from __future__ import annotations

from typing import TYPE_CHECKING, overload

if TYPE_CHECKING:
    from collections.abc import Iterable, Iterator, Mapping, Set
    from typing import Any, Final, Literal, NoReturn

    from docutils import nodes
    from typing_extensions import Self

    from sphinx.domains import Domain
    from sphinx.domains.c import CDomain
    from sphinx.domains.changeset import ChangeSetDomain
    from sphinx.domains.citation import CitationDomain
    from sphinx.domains.cpp import CPPDomain
    from sphinx.domains.index import IndexDomain
    from sphinx.domains.javascript import JavaScriptDomain
    from sphinx.domains.math import MathDomain
    from sphinx.domains.python import PythonDomain
    from sphinx.domains.rst import ReSTDomain
    from sphinx.domains.std import StandardDomain
    from sphinx.environment import BuildEnvironment
    from sphinx.ext.duration import DurationDomain
    from sphinx.ext.todo import TodoDomain


class _DomainsContainer:
    """Container for domain instances.

    This class is private, including its name, constructor, and all methods.
    Any or all of these will change without notice or warning in any release.

    The public interface is restricted to:

    * the``domains.['<domain-name>']`` mapping interface
    * the ``domains.<core-domain-name>`` attributes for core domains.
    * the `.get()``, ``.keys()``, ``.items()``, and ``.values()`` methods.

    Additionally, this class supports ``iter`` and ``len``,
    and provides membership testing via the ``in`` operator.
    """

    __slots__ = (
        '_domain_instances',
        'c_domain',
        'changeset_domain',
        'citation_domain',
        'cpp_domain',
        'index_domain',
        'javascript_domain',
        'math_domain',
        'python_domain',
        'restructuredtext_domain',
        'standard_domain',
    )

    #: First-party domains in :mod:`sphinx.domains`
    _core_domains: Final = frozenset({
        'std',
        # Language-specific domains
        'c',
        'cpp',
        'js',
        'py',
        'rst',
        # Other core domains
        'changeset',
        'citation',
        'index',
        'math',
    })

    @classmethod
    def _from_environment(cls, env: BuildEnvironment, /) -> Self:
        create_domains = env.app.registry.create_domains
        # Initialise domains
        if domains := {domain.name: domain for domain in create_domains(env)}:
            return cls(**domains)  # type: ignore[arg-type]

        return cls._from_environment_default(env=env)

    @classmethod
    def _from_environment_default(cls, *, env: BuildEnvironment) -> Self:
        """Return a default instance with every domain we require."""
        from sphinx.domains.c import CDomain
        from sphinx.domains.changeset import ChangeSetDomain
        from sphinx.domains.citation import CitationDomain
        from sphinx.domains.cpp import CPPDomain
        from sphinx.domains.index import IndexDomain
        from sphinx.domains.javascript import JavaScriptDomain
        from sphinx.domains.math import MathDomain
        from sphinx.domains.python import PythonDomain
        from sphinx.domains.rst import ReSTDomain
        from sphinx.domains.std import StandardDomain

        return cls(
            c=CDomain(env),
            changeset=ChangeSetDomain(env),
            citation=CitationDomain(env),
            cpp=CPPDomain(env),
            index=IndexDomain(env),
            js=JavaScriptDomain(env),
            math=MathDomain(env),
            py=PythonDomain(env),
            rst=ReSTDomain(env),
            std=StandardDomain(env),
        )

    def __init__(
        self,
        *,
        c: CDomain,
        cpp: CPPDomain,
        js: JavaScriptDomain,
        py: PythonDomain,
        rst: ReSTDomain,
        std: StandardDomain,
        changeset: ChangeSetDomain,
        citation: CitationDomain,
        index: IndexDomain,
        math: MathDomain,
        **domains: Domain,
    ) -> None:
        # All domains, including core.
        # Implemented as a dict for backwards compatibility.
        self._domain_instances: Mapping[str, Domain] = {
            'c': c,
            'changeset': changeset,
            'citation': citation,
            'cpp': cpp,
            'index': index,
            'js': js,
            'math': math,
            'py': py,
            'rst': rst,
            'std': std,
            **domains,
        }

        # Provide typed attributes for the core domains
        self.standard_domain: StandardDomain = std
        self.c_domain: CDomain = c
        self.cpp_domain: CPPDomain = cpp
        self.javascript_domain: JavaScriptDomain = js
        self.python_domain: PythonDomain = py
        self.restructuredtext_domain: ReSTDomain = rst
        self.changeset_domain: ChangeSetDomain = changeset
        self.citation_domain: CitationDomain = citation
        self.index_domain: IndexDomain = index
        self.math_domain: MathDomain = math

        for domain_name, domain in self._domain_instances.items():
            # invariant from ``_DomainsContainer._from_environment``
            if domain_name != domain.name:
                msg = f'Domain name mismatch in {domain!r}: {domain_name!r} != {domain.name!r}'
                raise ValueError(msg)

    def _setup(self) -> None:
        for domain in self._domain_instances.values():
            domain.setup()

    def _process_doc(
        self, env: BuildEnvironment, docname: str, document: nodes.document
    ) -> None:
        for domain in self._domain_instances.values():
            domain.process_doc(env, docname, document)

    def _clear_doc(self, docname: str) -> None:
        for domain in self._domain_instances.values():
            domain.clear_doc(docname)

    def _merge_domain_data(
        self, docnames: Set[str], domain_data: dict[str, Any]
    ) -> None:
        for domain_name, domain in self._domain_instances.items():
            domain.merge_domaindata(docnames, domain_data[domain_name])

    def _check_consistency(self) -> None:
        for domain in self._domain_instances.values():
            domain.check_consistency()

    def __contains__(self, key: str) -> bool:
        return key in self._domain_instances

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, _DomainsContainer):
            return NotImplemented
        return self._domain_instances == other._domain_instances

    def __setattr__(self, key: str, value: object) -> None:
        if key in self._core_domains:
            msg = f'{self.__class__.__name__!r} object does not support assignment to {key!r}'
            raise TypeError(msg)
        super().__setattr__(key, value)

    def __delattr__(self, key: str) -> None:
        if key in self._core_domains:
            msg = f'{self.__class__.__name__!r} object does not support deletion of {key!r}'
            raise TypeError(msg)
        super().__delattr__(key)

    # Mapping interface: builtin domains

    @overload
    def __getitem__(self, key: Literal['c']) -> CDomain: ...  # NoQA: E704

    @overload
    def __getitem__(self, key: Literal['cpp']) -> CPPDomain: ...  # NoQA: E704

    @overload
    def __getitem__(self, key: Literal['changeset']) -> ChangeSetDomain: ...  # NoQA: E704

    @overload
    def __getitem__(self, key: Literal['citation']) -> CitationDomain: ...  # NoQA: E704

    @overload
    def __getitem__(self, key: Literal['index']) -> IndexDomain: ...  # NoQA: E704

    @overload
    def __getitem__(self, key: Literal['js']) -> JavaScriptDomain: ...  # NoQA: E704

    @overload
    def __getitem__(self, key: Literal['math']) -> MathDomain: ...  # NoQA: E704

    @overload
    def __getitem__(self, key: Literal['py']) -> PythonDomain: ...  # NoQA: E704

    @overload
    def __getitem__(self, key: Literal['rst']) -> ReSTDomain: ...  # NoQA: E704

    @overload
    def __getitem__(self, key: Literal['std']) -> StandardDomain: ...  # NoQA: E704

    # Mapping interface: first-party domains

    @overload
    def __getitem__(self, key: Literal['duration']) -> DurationDomain: ...  # NoQA: E704

    @overload
    def __getitem__(self, key: Literal['todo']) -> TodoDomain: ...  # NoQA: E704

    # Mapping interface: third-party domains

    @overload
    def __getitem__(self, key: str) -> Domain: ...  # NoQA: E704

    def __getitem__(self, key: str) -> Domain:
        if domain := getattr(self, key, None):
            return domain
        return self._domain_instances[key]

    def __setitem__(self, key: str, value: Domain) -> NoReturn:
        msg = f'{self.__class__.__name__!r} object does not support item assignment'
        raise TypeError(msg)

    def __delitem__(self, key: str) -> NoReturn:
        msg = f'{self.__class__.__name__!r} object does not support item deletion'
        raise TypeError(msg)

    def __iter__(self) -> Iterator[str]:
        return iter(self._domain_instances.keys())

    def __len__(self) -> int:
        return len(self._domain_instances)

    def get(self, key: str, default: Domain | None = None) -> Domain | None:
        return self._domain_instances.get(key, default)

    def keys(self) -> Iterable[str]:
        return self._domain_instances.keys()

    def items(self) -> Iterable[tuple[str, Domain]]:
        return self._domain_instances.items()

    def values(self) -> Iterable[Domain]:
        return self._domain_instances.values()

    def sorted(self) -> Iterable[Domain]:
        for _domain_name, domain in sorted(self._domain_instances.items()):
            yield domain
