from __future__ import annotations

import sys
from typing import TYPE_CHECKING

if sys.version_info >= (3, 10):
    from typing import TypeAlias
else:
    from typing_extensions import TypeAlias

if sys.version_info >= (3, 11):
    from typing import Self
elif TYPE_CHECKING:
    from typing_extensions import Self
else:
    Self = object


__all__ = ["Self", "TypeAlias"]


def __dir__() -> list[str]:
    return __all__
