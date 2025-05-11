from pathlib import Path
import time

# Test path to >1M reads BAM file.
path_to_bam = Path(
    r"C:\Users\nogtr\data\basecalled_seqdata\250324_ecoli_WT_RCCmix1_v17D1to8\calls.bam"
)
print(path_to_bam)
print(path_to_bam.exists())
print(path_to_bam.is_file())

import bampy as bp

f = bp.BamReader(str(path_to_bam))

start_time = time.time()
print("start_time:", start_time)
count = 0
for record in f:
    count += 1
    if count % 100000 == 0:
        print(count)

elapsed_time = time.time() - start_time
print("elapsed_time:", elapsed_time)
