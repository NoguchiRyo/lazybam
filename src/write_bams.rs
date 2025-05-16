use noodles::bam;
use noodles::sam;
use sam::alignment::io::Write;
use std::{fs::File, path::Path};

/// Write a chunk of alignments to a BAM file.
///
/// A convenience wrapper that **optionally coordinate‑sorts** a mutable slice of
/// `RecordBuf`s and streams them to disk in BGZF‑compressed BAM format.
///
/// Parameters
/// ----------
/// * `header` – Parsed SAM header describing the reference sequence dictionary
///   and any metadata that must precede the records in a valid BAM.
/// * `records` – Mutable vector of owned alignment records (`RecordBuf`).  
///   When `sort == true` the vector is reordered in‑place; otherwise record
///   order is preserved.
/// * `path` – Destination file path. The file is created (clobbering any
///   existing file) and wrapped in a `noodles_bam::io::Writer`.
/// * `sort` – If `true`, records are coordinate‑sorted
///   (`reference_sequence_id`, then `alignment_start`) **before** writing.  
///   Sorting is done with `rayon::par_sort_unstable_by`, exploiting all CPU
///   cores.
///
/// Returns
/// -------
/// * `std::io::Result<()>` – `Ok(())` on success or any I/O/encoding error
///   propagated from `noodles` internals.
///
/// Notes
/// -----
/// * Uses `Writer::new` → BGZF compression (samtools‑compatible).
/// * Explicitly calls `writer.try_finish()` to flush BGZF blocks before the
///   `File` is dropped, guaranteeing on‑disk consistency even if the caller
///   immediately re‑opens the file.
///
/// Example
/// -------
/// ```
/// write_chunk(&header, &mut chunk, "chunk_001.bam", true)?;
/// ```
pub fn write_chunk<P>(
    header: &sam::Header,
    records: &mut Vec<sam::alignment::RecordBuf>,
    path: P,
    sort: bool,
) -> std::io::Result<()>
where
    P: AsRef<Path>,
{
    // ── 1. Optional parallel coordinate sort ────────────────────────────────
    if sort && records.len() > 1 {
        records.sort_unstable_by(|a, b| coord_cmp(a, b));
    }

    // ── 2. Initialise BGZF‑compressed BAM writer ────────────────────────────
    let file = File::create(path)?;
    let mut writer = bam::io::Writer::new(file);

    // ── 3. Emit header and alignment records ────────────────────────────────
    writer.write_header(header)?;
    for rec in records {
        writer.write_alignment_record(header, rec)?;
    }

    // ── 4. Flush BGZF blocks to disk ────────────────────────────────────────
    writer.try_finish()?;
    Ok(())
}

/// Total ordering for coordinate sort used by `write_chunk`.
///
/// Ordering rules
/// --------------
/// 1. **Mapped vs. unmapped** – mapped reads (`Some(id)`) precede unmapped
///    reads (`None`).
/// 2. **Reference sequence ID** – ascending numeric `ref_id`.
/// 3. **Alignment start** – ascending 1‑based coordinate within the reference.
fn coord_cmp(a: &sam::alignment::RecordBuf, b: &sam::alignment::RecordBuf) -> std::cmp::Ordering {
    use std::cmp::Ordering::*;
    match (a.reference_sequence_id(), b.reference_sequence_id()) {
        // Both mapped: compare reference, then position
        (Some(ra), Some(rb)) => match ra.cmp(&rb) {
            Equal => a
                .alignment_start()
                .expect("Invalid alignment start")
                .cmp(&b.alignment_start().expect("Invalid alignment start")),
            non_eq => non_eq,
        },
        // Mapped < unmapped
        (Some(_), None) => Less,
        (None, Some(_)) => Greater,
        // Both unmapped
        (None, None) => Equal,
    }
}
