from .lazybam import (
    BamReader,
    write_chunk_py,
    merge_chunks_py,
    PyBamRecord,
    RecordOverride,
    PyRecordBuf,
    write_recordbuf_chunk_py,
)
from .header import BamHeader

__all__ = [
    "BamReader",
    "write_chunk_py",
    "merge_chunks_py",
    "BamHeader",
    "PyBamRecord",
    "RecordOverride",
    "PyRecordBuf",
    "write_recordbuf_chunk_py",
]


def _get_header(self: BamReader):
    # self._header は RustBamReader 側で定義された bytes 型の属性
    return BamHeader.from_bytes(self._header)


# property としてクラスに追加
BamReader.header = property(_get_header)  #  type: ignore

__doc__ = lazybam.__doc__
