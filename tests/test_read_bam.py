from pathlib import Path

path_to_bam = Path(__file__).parent / "data" / "test_reads.bam"
print(path_to_bam)
print(path_to_bam.exists())
print(path_to_bam.is_file())

import bampy as bp

f = bp.BamReader(str(path_to_bam))

for record in f:
    print(record)
    print(record.qname)
    print(record.seq)
    print(record.qual)
    print(record.qual.shape)
    print(record.tags)
    print(record.cigar)
