from pathlib import Path

path_to_bam = Path(__file__).parent / "data" / "test_reads.bam"
print(path_to_bam)
print(path_to_bam.exists())
print(path_to_bam.is_file())

import bampy._bampy as bp

f = bp.BamReader(str(path_to_bam), chunk_size=1000)

print(f.header.decode("utf-8"))

for records in f:
    for record in records:
        print(record)
        print(record.qname)
        print(record.seq)
        print(record.qual)
        print(record.qual.shape)
        print(record.tags)
        print(record.cigar)
