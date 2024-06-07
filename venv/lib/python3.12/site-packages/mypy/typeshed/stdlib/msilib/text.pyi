import sys

if sys.platform == "win32":
    ActionText: list[tuple[str, str, str | None]]
    UIText: list[tuple[str, str | None]]
    dirname: str
    tables: list[str]
