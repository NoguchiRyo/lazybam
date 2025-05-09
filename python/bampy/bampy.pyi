# bampy.pyi  ------------------------------------------------------------
from __future__ import annotations

from typing import Any, Iterator, List, Tuple

import numpy as np
import numpy.typing as npt

__all__: list[str]  # ← 好みで

# ----------------------------------------------------------------------
class PyBamRecord:
    # ── public attributes ------------------------------------------------
    qname: str
    flag: int
    pos: int
    len: int  # template length
    mapq: int
    tags: List[Tuple[str, Any]]  # 各タグは数値 / str / NumPy 配列など混在

    # ── getters (読み取り専用プロパティ) ----------------------------------
    @property
    def seq(self) -> str: ...
    @property
    def qual(self) -> npt.NDArray[np.uint8]: ...
    @property
    def cigar(self) -> npt.NDArray[np.uint8]: ...

# ----------------------------------------------------------------------
class BamReader(Iterator[PyBamRecord]):
    def __init__(self, path: str) -> None: ...

    # ── context‑manager --------------------------------------------------
    def __enter__(self) -> BamReader: ...
    def __exit__(
        self,
        exc_type: Any,
        exc_val: Any,
        traceback: Any,
    ) -> None: ...

    # ── iterator ---------------------------------------------------------
    def __iter__(self) -> BamReader: ...
    def __next__(self) -> PyBamRecord: ...

    # ── other properties -------------------------------------------------
    @property
    def header(self) -> bytes: ...

# ----------------------------------------------------------------------
def sum_as_string(a: int, b: int) -> str: ...
