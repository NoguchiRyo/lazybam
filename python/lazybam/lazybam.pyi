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
        qname: Optional[str] = None,
        seq: Optional[str] = None,
        qual: Optional[List[int]] = None,
        reference_sequence_id: Optional[int] = None,
        cigar: Optional[List[Tuple[int, int]]] = None,
        alignment_start: Optional[int] = None,
        tags: Optional[List[Tuple[str, Any]]] = None,
        mapping_quality: Optional[int] = None,
    ) -> None: ...
    @property
    def qname(self) -> Optional[str]: ...
    @qname.setter
    def qname(self, name: str) -> None: ...
    @property
    def seq(self) -> Optional[str]: ...
    @seq.setter
    def seq(self, sequence: str) -> None: ...
    @property
    def qual(self) -> Optional[List[int]]: ...
    @qual.setter
    def qual(self, quality: List[int]) -> None: ...
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
    @property
    def alignment_start(self) -> Optional[int]: ...
    @alignment_start.setter
    def alignment_start(self, pos: int) -> None: ...
    @property
    def mapping_quality(self) -> Optional[int]: ...
    @mapping_quality.setter
    def mapping_quality(self, mapq: int) -> None: ...

class PyBamRecord:
    # ── public attributes ------------------------------------------------
    qname: str
    flag: int
    pos: int
    len: int  # template length
    mapq: int
    rid: int

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
    def get_field_by_tag(self, tag: str) -> Any: ...

class PyRecordBuf:
    def __init__(
        self,
        qname: str,
        seq: str,
        qual: List[int],
        reference_sequence_id: Optional[int] = None,
        cigar: Optional[List[Tuple[int, int]]] = None,
        alignment_start: Optional[int] = None,
        mapping_quality: Optional[int] = None,
        tags: Optional[List[Tuple[str, Any]]] = None,
    ) -> None: ...

class BamReader:
    def __init__(
        self, path: str, chunk_size: Optional[int] = None, region: Optional[str] = None
    ) -> None: ...

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
def write_recordbuf_chunk_py(
    header_bytes: bytes,
    records: List[PyRecordBuf],
    out_bam: str,
    sort: bool,
) -> None: ...
