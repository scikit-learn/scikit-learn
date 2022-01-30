from typing import Any

from ..cmd import Command

_Reporter = Any  # really docutils.utils.Reporter

# Only defined if docutils is installed.
class SilentReporter(_Reporter):  # type: ignore
    messages: Any
    def __init__(
        self,
        source,
        report_level,
        halt_level,
        stream: Any | None = ...,
        debug: int = ...,
        encoding: str = ...,
        error_handler: str = ...,
    ) -> None: ...
    def system_message(self, level, message, *children, **kwargs): ...

HAS_DOCUTILS: bool

class check(Command):
    description: str
    user_options: Any
    boolean_options: Any
    restructuredtext: int
    metadata: int
    strict: int
    def initialize_options(self) -> None: ...
    def finalize_options(self) -> None: ...
    def warn(self, msg): ...
    def run(self) -> None: ...
    def check_metadata(self) -> None: ...
    def check_restructuredtext(self) -> None: ...
