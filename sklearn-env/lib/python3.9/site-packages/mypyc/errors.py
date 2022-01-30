from typing import List

import mypy.errors


class Errors:
    def __init__(self) -> None:
        self.num_errors = 0
        self.num_warnings = 0
        self._errors = mypy.errors.Errors()

    def error(self, msg: str, path: str, line: int) -> None:
        self._errors.report(line, None, msg, severity='error', file=path)
        self.num_errors += 1

    def note(self, msg: str, path: str, line: int) -> None:
        self._errors.report(line, None, msg, severity='note', file=path)

    def warning(self, msg: str, path: str, line: int) -> None:
        self._errors.report(line, None, msg, severity='warning', file=path)
        self.num_warnings += 1

    def new_messages(self) -> List[str]:
        return self._errors.new_messages()

    def flush_errors(self) -> None:
        for error in self.new_messages():
            print(error)
