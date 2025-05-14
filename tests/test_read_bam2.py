from pathlib import Path
import time

# Test path to >1M reads BAM file.
path_to_bam = Path(
    r"C:\Users\nogtr\data\basecalled_seqdata\250324_ecoli_WT_RCCmix1_v17D1to8\calls.bam"
)
print(path_to_bam)
print(path_to_bam.exists())
print(path_to_bam.is_file())

import lazybam as lb

f = lb.BamReader(str(path_to_bam), chunk_size=1)

start_time = time.time()
print("start_time:", start_time)
count = 0
for records in f:
    for record in records:
        # print(record)
        # print(record.qname)
        # print(record.seq)
        # print(record.qual)
        # print(record.qual.shape)
        # print(record.tags)
        # print(record.cigar)
        count += 1
        if count % 100000 == 0:
            print(count)

elapsed_time = time.time() - start_time
print("elapsed_time:", elapsed_time)
