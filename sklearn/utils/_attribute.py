# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

from dataclasses import dataclass, field


@dataclass
class FittedAttribute:
    """Descriptor to raise a better error when attribute is not defined."""

    condition: str = ""
    name: str = field(init=False)

    def __get__(self, instance, owner):
        if self.condition:
            msg = f"Call `fit` with {self.condition} to define '{self.name}'"
        else:
            msg = f"Call `fit` {self.name} to define '{self.name}'"
        raise AttributeError(msg)

    def __set_name__(self, owner, name):
        self.name = name
