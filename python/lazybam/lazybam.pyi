from __future__ import annotations

from typing import Any, List, Optional, Tuple

import numpy as np  # type: ignore

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from .header import BamHeader

__doc__: str = "Rust powered BAM reader built on noodles + PyO3"

class PyKind:
    Match: PyKind
    Insertion: PyKind
    Deletion: PyKind
    Skip: PyKind
    SoftClip: PyKind
    HardClip: PyKind
    Pad: PyKind
    SequenceMatch: PyKind
    SequenceMismatch: PyKind

class RecordOverride:
    def __init__(
        self,
        reference_sequence_id: Optional[int] = None,
        cigar: Optional[List[Tuple[int, int]]] = None,
        tags: Optional[List[Tuple[str, Any]]] = None,
    ) -> None: ...
    @property
    def reference_sequence_id(self) -> Optional[int]: ...
    @reference_sequence_id.setter
    def reference_sequence_id(self, rid: int) -> None: ...
    @property
    def cigar(self) -> Optional[List[Tuple[int, int]]]: ...
    @cigar.setter
    def cigar(self, cigar_list: List[Tuple[int, int]]) -> None: ...
    @property
    def tags(self) -> List[Tuple[str, Any]]: ...
    @tags.setter
    def tags(self, vals: List[Tuple[str, Any]]) -> None: ...

class PyBamRecord:
    # ── public attributes ------------------------------------------------
    qname: str
    flag: int
    pos: int
    len: int  # template length
    mapq: int

    # ── getters (read-only properties) ----------------------------------
    @property
    def seq(self) -> str: ...
    @property
    def qual(self) -> List[int]: ...
    @property
    def cigar(self) -> List[Tuple[int, int]]: ...
    @property
    def tags(self) -> List[Tuple[str, Any]]: ...
    def set_record_override(self, record_override: RecordOverride) -> None: ...

class BamReader:
    def __init__(self, path: str, chunk_size: Optional[int] = None) -> None: ...

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
    def __next__(self) -> List[PyBamRecord]: ...

    # ── other properties -------------------------------------------------
    @property
    def _header(self) -> bytes: ...
    @property
    def header(self) -> BamHeader: ...

# Writing functions
def write_chunk_py(
    header_bytes: bytes,
    records: List[PyBamRecord],
    out_bam: str,
    sort: bool,
) -> None: ...
def merge_chunks_py(
    header_bytes: bytes,
    chunks: List[str],
    out_bam: str,
    sort: bool,
) -> None: ...
