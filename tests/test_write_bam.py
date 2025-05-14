from pathlib import Path

path_to_bam = Path(__file__).parent / "data" / "test_reads.bam"
print(path_to_bam)
print(path_to_bam.exists())
print(path_to_bam.is_file())

import lazybam as lb

f = lb.BamReader(str(path_to_bam), chunk_size=1000)

print(f.header.decode("utf-8"))
record_list: list[lb.PyBamRecord] = []
for records in f:
    for record in records:
        record_list.append(record)

print("record length:", len(record_list))

chunk_path = path_to_bam.parent / "test_reads_out_chunk.bam"
out_path = path_to_bam.parent / "test_reads_out.bam"
lb.write_chunk_py(
    header_bytes=f.header,
    records=record_list,
    out_bam=str(chunk_path),
    sort=True,
)
lb.merge_chunks_py(
    header_bytes=f.header,
    chunks=[str(chunk_path)],
    out_bam=str(out_path),
    sort=True,
)
