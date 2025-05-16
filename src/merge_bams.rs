//! Utilities for merging multiple chunk BAM files into a single BAM and,
//! if requested, creating the corresponding coordinate–sorted BAI index.
//!
//! The public entry point is [`merge_chunks`]. All helpers are kept private
//! to minimize surface area.

use noodles::csi::binning_index;
use noodles::{bam, bgzf, sam};
use sam::alignment::{io::Write as _, RecordBuf};
use std::{
    cmp::Ordering,
    collections::BinaryHeap,
    fs::File,
    path::{Path, PathBuf},
};

/// Merge a set of temporary chunk BAMs into `out_bam`.
///
/// When `sort == true` the routine performs a *global* k‑way merge that
/// preserves coordinate order **and** writes a `*.bai` index side‑by‑side.
/// Otherwise the chunks are concatenated in the order given and no index
/// is produced.
///
/// ### Parameters
/// * **`header`** – A [`sam::Header`] identical for all chunks and the final BAM.
/// * **`chunks`** – Paths to the temporary chunk files, typically produced by
///   your earlier `write_chunk`.
/// * **`out_bam`** – Destination path for the merged BAM (existing file is
///   overwritten).
/// * **`sort`** – `true` ⇒ k‑way coordinate merge + BAI generation.  
///                `false` ⇒ raw concatenation, no index.
///
/// ### Errors
/// Propagates any I/O or codec error returned by `noodles` crates.
pub fn merge_chunks(
    header: &sam::Header,
    chunks: &[PathBuf],
    out_bam: &Path,
    sort: bool,
) -> std::io::Result<()> {
    // println!(
    //     "merge_chunks: starting with {} chunks, sort={}...",
    //     chunks.len(),
    //     sort
    // );

    // ── 1. Set up the output writer ───────────────────────────────────────
    let file = File::create(out_bam)?;
    let mut writer = bam::io::Writer::new(file);
    // println!("merge_chunks: writing header to {:?}", out_bam);
    writer.write_header(header)?;

    // ── 2. Fast path: just concatenate chunks, no global sort ─────────────
    if !sort {
        for p in chunks {
            // println!("merge_chunks: concatenating chunk {:?}...", p);
            copy_records(header, p, &mut writer)?;
        }
        // println!("merge_chunks: flushing writer (concatenate path)");
        writer.try_finish()?; // flush BGZF
        return Ok(());
    }

    // ── 3. Build k‑way merge structures ───────────────────────────────────
    // println!("merge_chunks: opening chunk readers...");
    let mut readers: Vec<bam::io::Reader<bgzf::Reader<File>>> = chunks
        .iter()
        .map(|p| {
            // println!("merge_chunks: opening chunk {:?}", p);
            bam::io::reader::Builder::default()
                .build_from_path(p)
                .and_then(|mut r| {
                    r.read_header()?; // skip chunk header
                    Ok(r)
                })
        })
        .collect::<Result<_, _>>()?;

    /// Lexicographic key used to order records in the heap.
    fn key(rec: &RecordBuf) -> (Option<usize>, usize, bool) {
        (
            rec.reference_sequence_id(),
            rec.alignment_start()
                .unwrap_or_else(|| noodles::core::Position::try_from(0).unwrap())
                .get(),
            rec.reference_sequence_id().is_none(),
        )
    }

    /// Heap entry: the current record for a given chunk.
    struct Entry {
        chunk_idx: usize,
        rec: RecordBuf,
    }

    // Implement ordering so that BinaryHeap becomes a min‑heap by coordinate.
    impl Eq for Entry {}
    impl PartialEq for Entry {
        fn eq(&self, other: &Self) -> bool {
            key(&self.rec) == key(&other.rec)
        }
    }
    impl Ord for Entry {
        fn cmp(&self, other: &Self) -> Ordering {
            // Invert comparison to make BinaryHeap act as a min-heap by coordinate
            key(&other.rec).cmp(&key(&self.rec))
        }
    }
    impl PartialOrd for Entry {
        fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
            Some(self.cmp(other))
        }
    }

    let mut heap: BinaryHeap<Entry> = BinaryHeap::new();
    // println!("merge_chunks: priming heap with one record per chunk...");

    // Prime the heap with one record from each reader.
    for (i, rdr) in readers.iter_mut().enumerate() {
        let mut rec = RecordBuf::default();
        if rdr.read_record_buf(header, &mut rec)? != 0 {
            // println!(
            //     "merge_chunks: chunk {} yields first record ref_id={:?}, start={:?}",
            //     i,
            //     rec.reference_sequence_id(),
            //     rec.alignment_start()
            // );
            heap.push(Entry { chunk_idx: i, rec });
        } else {
            // println!("merge_chunks: chunk {} was empty", i);
        }
    }

    // ── 4. Perform streaming k‑way merge ──────────────────────────────────
    // println!("merge_chunks: entering merge loop...");
    // let mut count = 0usize;
    while let Some(Entry { chunk_idx, rec }) = heap.pop() {
        writer.write_alignment_record(header, &rec)?;
        // count += 1;
        // if count % 100_000 == 0 {
        //     println!("merge_chunks: wrote {} records...", count);
        // }

        let rdr = &mut readers[chunk_idx];
        let mut next = RecordBuf::default();
        if rdr.read_record_buf(header, &mut next)? != 0 {
            heap.push(Entry {
                chunk_idx,
                rec: next,
            });
        }
    }
    // println!("merge_chunks: merge loop done, total {} records", count);

    writer.try_finish()?; // flush BGZF blocks
                          // println!("merge_chunks: writer flushed");

    // ── 5. Build BAI index in a second pass over the merged BAM ───────────
    if let Some(idx_path) = out_bam.to_str().map(|s| format!("{s}.bai")) {
        // println!("merge_chunks: building BAI index at {:?}", idx_path);
        build_bai_index(header, out_bam, Path::new(&idx_path))?;
        // println!("merge_chunks: BAI index built");
    }

    // println!("merge_chunks: done");
    Ok(())
}

// ==========================================================================
// Helper functions
// ==========================================================================

/// Stream‑copy all records from `src_bam` into `writer`.
///
/// Used when `sort == false` to avoid the k‑way merge overhead.
fn copy_records(
    header: &sam::Header,
    src_bam: &Path,
    writer: &mut bam::io::Writer<bgzf::Writer<File>>,
) -> std::io::Result<()> {
    println!("copy_records: open {:?}", src_bam);
    let mut reader = bam::io::reader::Builder::default().build_from_path(src_bam)?;
    reader.read_header()?; // discard chunk header

    let mut rec = RecordBuf::default();
    let mut cnt = 0;
    while reader.read_record_buf(header, &mut rec)? != 0 {
        writer.write_alignment_record(header, &rec)?;
        cnt += 1;
        if cnt % 100_000 == 0 {
            println!("copy_records: wrote {} records from {:?}", cnt, src_bam);
        }
    }
    println!("copy_records: done {:?}, total {} records", src_bam, cnt);
    Ok(())
}

/// Build a BAI index by scanning `bam_path` once and serialising to `bai_path`.
///
/// Uses the low‑level `binning_index::Indexer` so it works across all released
/// versions of `noodles`.
fn build_bai_index(header: &sam::Header, bam_path: &Path, bai_path: &Path) -> std::io::Result<()> {
    println!("build_bai_index: scanning {:?}", bam_path);
    use binning_index::index::reference_sequence::bin::Chunk;

    let mut reader = bam::io::reader::Builder::default().build_from_path(bam_path)?;
    reader.read_header()?; // header equals `header`

    // Default indexer implements the BAM‑specific 14‑level / 5‑granularity binning.
    let mut indexer = binning_index::Indexer::<
        binning_index::index::reference_sequence::index::LinearIndex,
    >::default();

    let mut record = RecordBuf::default();
    let mut chunk_start = reader.get_ref().virtual_position();
    let mut idx_cnt = 0;

    while reader.read_record_buf(header, &mut record)? != 0 {
        let chunk_end = reader.get_ref().virtual_position();

        let alignment_ctx = match (
            record.reference_sequence_id(),
            record.alignment_start(),
            record.alignment_end(),
        ) {
            (Some(rid), Some(start_pos), Some(end_pos)) => {
                let mapped = !record
                    .flags()
                    .contains(sam::alignment::record::Flags::UNMAPPED);
                // println!(
                //     "build_bai_index: add_record rid={} start={:?} end={:?}",
                //     rid, start_pos, end_pos
                // );
                Some((rid, start_pos, end_pos, mapped))
            }
            _ => None, // unplaced or unmapped
        };

        indexer.add_record(alignment_ctx, Chunk::new(chunk_start, chunk_end))?;
        idx_cnt += 1;
        if idx_cnt % 100_000 == 0 {
            println!("build_bai_index: processed {} records", idx_cnt);
        }

        chunk_start = reader.get_ref().virtual_position();
    }

    println!(
        "build_bai_index: finalizing index after {} records",
        idx_cnt
    );
    let index = indexer.build(header.reference_sequences().len());

    println!("build_bai_index: writing index to {:?}", bai_path);
    let bai_file = File::create(bai_path)?;
    let mut bai_writer = bam::bai::io::Writer::new(bai_file);
    bai_writer.write_index(&index)?;
    println!("build_bai_index: done");

    Ok(())
}
