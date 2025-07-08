# Lazybam

## Introduction

A fast, memory-efficient, and OS-independent Python library for reading and writing BAM files. Built on Rust's noodles library with PyO3 bindings, lazybam provides high-performance BAM file processing with a simple Python API across Windows, macOS, and Linux.

## Installation

```bash
pip install lazybam
```

## Usage

### Basic BAM File Reading

```python
import lazybam as lb

# Read a BAM file with chunked processing
with lb.BamReader("path/to/file.bam", chunk_size=1000) as reader:
    # Access header information
    header = reader.header
    print(header)
    
    # Iterate through records in chunks
    for chunk in reader:
        for record in chunk:
            print(f"Read name: {record.qname}")
            print(f"Sequence: {record.seq}")
            print(f"Quality: {record.qual}")
            print(f"CIGAR: {record.cigar}")
            print(f"Tags: {record.tags}")
```

### Region-Specific Queries

```python
# Read specific genomic regions
reader = lb.BamReader("path/to/file.bam", region="chr1:1000-2000")
for chunk in reader:
    for record in chunk:
        ref_name = reader.header.get_ref_name(record.rid)
        print(f"Reference: {ref_name}, Position: {record.pos}")
```

### Working with Headers

```python
# Access and modify BAM headers
header = reader.header
print(f"Reference name for RID 0: {header.get_ref_name(0)}")

# Add program information
header.add_program(ID="myprogram", VN="1.0", CL="mycommand")

# Change sorting order
header.change_SO_tag("coordinate")
```

### Writing BAM Files

```python
# Collect records for writing
records = []
with lb.BamReader("input.bam") as reader:
    header_bytes = reader.header.to_bytes()
    for chunk in reader:
        records.extend(chunk)

# Write records to a new BAM file
lb.write_chunk_py(
    header_bytes=header_bytes,
    records=records,
    out_bam="output.bam",
    sort=True
)
```

### Creating Custom Records

```python
# Create records from scratch
record = lb.PyRecordBuf(
    qname="read1",
    seq="ATCGATCG",
    qual=[30, 30, 30, 30, 30, 30, 30, 30],
    reference_sequence_id=0,
    alignment_start=100,
    mapping_quality=60,
    cigar=[(0, 8)],  # 8M
    tags=[("NM", 0)]
)

# Write custom records
lb.write_recordbuf_chunk_py(
    header_bytes=header_bytes,
    records=[record],
    out_bam="custom.bam",
    sort=True
)
```

### Merging BAM Files

```python
# Merge multiple BAM chunks
chunk_files = ["chunk1.bam", "chunk2.bam", "chunk3.bam"]
lb.merge_chunks_py(
    header_bytes=header_bytes,
    chunks=chunk_files,
    out_bam="merged.bam",
    sort=True  # When True, also generates a .bai index file
)
```

### Modifying Records

```python
# Override specific record fields
override = lb.RecordOverride(
    qname="new_name",
    seq="AAAATTTTCCCCGGGG",
    mapping_quality=30
)

# Apply override to a record
record.set_record_override(override)
```
