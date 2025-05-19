from pathlib import Path

path_to_bam = Path(
    r"C:\Users\nogtr\data\sequence_alignment_data\250217_ecoli_ttcA_total_LB_BW_sta_v17Db1\250217_ecoli_ttcA_total_LB_BW_sta_v17Db1_sorted.bam"
)
print(path_to_bam)
print(path_to_bam.exists())
print(path_to_bam.is_file())

import lazybam as lb

f = lb.BamReader(str(path_to_bam), chunk_size=1, region="Arg4")
header = f.header

count = 0
for records in f:
    for record in records:
        # print(record.rid)
        print(header.get_ref_name(record.rid))
        count += 1
    if count > 100:
        break
