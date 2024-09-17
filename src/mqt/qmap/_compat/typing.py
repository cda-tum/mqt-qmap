# noqa: A005
from __future__ import annotations

import sys
from typing import TYPE_CHECKING

if sys.version_info < (3, 10):
    from typing_extensions import TypeAlias
else:
    from typing import TypeAlias

if sys.version_info < (3, 11):
    if TYPE_CHECKING:
        from typing_extensions import Self
    else:
        Self = object
else:
    from typing import Self

__all__ = ["Self", "TypeAlias"]


def __dir__() -> list[str]:
    return __all__
