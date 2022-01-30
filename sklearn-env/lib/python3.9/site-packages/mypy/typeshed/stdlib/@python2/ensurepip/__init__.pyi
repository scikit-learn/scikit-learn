def version() -> str: ...
def bootstrap(
    root: str | None = ...,
    upgrade: bool = ...,
    user: bool = ...,
    altinstall: bool = ...,
    default_pip: bool = ...,
    verbosity: int = ...,
) -> None: ...
