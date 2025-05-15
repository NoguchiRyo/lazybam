use numpy::PyArray1;
use pyo3::prelude::*;
use pyo3::IntoPyObjectExt;

use noodles::sam::alignment::record::data::field::value::Array;
use noodles::sam::alignment::record_buf::{Cigar, Data, QualityScores, Sequence as SeqBuf};
use noodles::sam::alignment::{
    record::{Flags, MappingQuality},
    RecordBuf,
};
use noodles::{bam, core::Position, sam};
use sam::alignment::record::cigar::op::Op;
use sam::alignment::record::data::field::Value as BamValue;
use sam::alignment::record::Cigar as _;

use crate::record_override::RecordOverride;

#[pyclass]
#[derive(Clone, Copy, Debug)]
pub enum PyKind {
    Match = 0,
    Insertion,
    Deletion,
    Skip,
    SoftClip,
    HardClip,
    Pad,
    SequenceMatch,
    SequenceMismatch,
}

impl From<sam::alignment::record::cigar::op::Kind> for PyKind {
    fn from(k: sam::alignment::record::cigar::op::Kind) -> Self {
        use sam::alignment::record::cigar::op::Kind::*;
        match k {
            Match => PyKind::Match,
            Insertion => PyKind::Insertion,
            Deletion => PyKind::Deletion,
            Skip => PyKind::Skip,
            SoftClip => PyKind::SoftClip,
            HardClip => PyKind::HardClip,
            Pad => PyKind::Pad,
            SequenceMatch => PyKind::SequenceMatch,
            SequenceMismatch => PyKind::SequenceMismatch,
        }
    }
}

#[pyclass]
pub struct PyBamRecord {
    record: bam::Record,
    record_override: Option<RecordOverride>,
}

impl PyBamRecord {
    pub fn from_record(record: bam::Record) -> Self {
        Self {
            record,
            record_override: None,
        }
    }

    /// Convert to RecordBuf, applying overrides
    pub fn to_record_buf(&self) -> anyhow::Result<RecordBuf> {
        // sequence & quality
        let seq = SeqBuf::from(self.seq().into_bytes());
        let qual = QualityScores::from(self.record.quality_scores().as_ref().to_vec());
        let mut data = Data::try_from(self.record.data()).unwrap_or_default();

        let position = match self.record.alignment_start() {
            Some(Ok(pos)) => pos,
            Some(Err(_)) => {
                return Err(anyhow::anyhow!("Invalid alignment start position"));
            }
            None => Position::try_from(0).unwrap(),
        };

        let mut cigar_vec: Vec<Op> = self.record.cigar().iter().filter_map(Result::ok).collect();
        let mut ref_id = match self.record.reference_sequence_id() {
            Some(Ok(rid)) => rid,
            Some(Err(_)) => {
                return Err(anyhow::anyhow!("Invalid reference sequence ID"));
            }
            None => 0,
        };
        if let Some(r_override) = &self.record_override {
            for (tag, value) in &r_override.tags {
                data.insert(*tag, value.clone());
            }
            if let Some(cigar) = &r_override.cigar {
                cigar_vec = cigar.iter().filter_map(Result::ok).collect();
            }

            if let Some(reference_id) = &r_override.reference_sequence_id {
                ref_id = *reference_id as usize;
            }
        }
        // builder
        let b = RecordBuf::builder()
            .set_name(self.qname())
            .set_flags(Flags::from_bits_truncate(self.flag()))
            .set_cigar(Cigar::from(cigar_vec))
            .set_sequence(seq)
            .set_quality_scores(qual)
            .set_data(data)
            .set_reference_sequence_id(ref_id)
            .set_alignment_start(position)
            .set_mapping_quality(MappingQuality::try_from(self.mapq()).unwrap());
        Ok(b.build())
    }
}

#[pymethods]
impl PyBamRecord {
    #[setter]
    fn set_record_override(&mut self, override_: RecordOverride) {
        self.record_override = Some(override_);
    }
    // ── getters and setters ────────────────────────────────────────────
    #[getter]
    fn qname(&self) -> String {
        self.record
            .name()
            .map(|b| b.to_string())
            .unwrap_or_default()
    }
    #[getter]
    fn flag(&self) -> u16 {
        u16::from(self.record.flags())
    }
    #[getter]
    fn pos(&self) -> i64 {
        self.record
            .alignment_start()
            .and_then(|r| r.ok())
            .map(|p| usize::from(p) as i64)
            .unwrap_or(-1)
    }
    #[getter]
    fn mapq(&self) -> u8 {
        self.record
            .mapping_quality()
            .map(|mq| u8::from(mq))
            .unwrap_or(255)
    }
    #[getter]
    fn len(&self) -> usize {
        self.record.template_length().abs() as usize
    }

    #[getter]
    fn seq(&self) -> String {
        self.record.sequence().iter().map(|b| b as char).collect()
    }
    #[getter]
    fn qual(&self) -> Vec<usize> {
        self.record
            .quality_scores()
            .as_ref()
            .iter()
            .map(|&b| b as usize)
            .collect()
    }

    #[getter]
    fn cigar(&self) -> Vec<(u32, u32)> {
        let ops: Vec<(u32, u32)> = self
            .record
            .cigar()
            .iter()
            .filter_map(Result::ok)
            .map(|op| (op.kind() as u32, op.len() as u32))
            .collect();
        return ops;
    }

    #[getter]
    fn tags<'py>(&self, py: Python<'py>) -> Vec<(String, PyObject)> {
        // override がなければ元の record.data() から構築
        let mut vec = Vec::new();
        for field in self.record.data().iter().filter_map(Result::ok) {
            let key = String::from_utf8_lossy(field.0.as_ref()).into_owned();
            let obj = match field.1 {
                BamValue::Int8(n) => (n as i32).into_py_any(py).unwrap(),
                BamValue::UInt8(n) => (n as u32).into_py_any(py).unwrap(),
                BamValue::Int16(n) => (n as i32).into_py_any(py).unwrap(),
                BamValue::UInt16(n) => (n as u32).into_py_any(py).unwrap(),
                BamValue::Int32(n) => (n as i32).into_py_any(py).unwrap(),
                BamValue::UInt32(n) => (n as u32).into_py_any(py).unwrap(),
                BamValue::Float(f) => (f as f64).into_py_any(py).unwrap(),
                BamValue::Character(c) => c.to_string().into_py_any(py).unwrap(),
                BamValue::String(bs) => String::from_utf8_lossy(bs)
                    .into_owned()
                    .into_py_any(py)
                    .unwrap(),
                BamValue::Array(arr) => match arr {
                    Array::UInt8(a) => {
                        PyArray1::from_vec(py, a.iter().filter_map(|r| r.ok()).collect())
                            .into_py_any(py)
                            .unwrap()
                    }
                    Array::Int8(a) => {
                        PyArray1::from_vec(py, a.iter().filter_map(|r| r.ok()).collect())
                            .into_py_any(py)
                            .unwrap()
                    }
                    Array::Int16(a) => {
                        PyArray1::from_vec(py, a.iter().filter_map(|r| r.ok()).collect())
                            .into_py_any(py)
                            .unwrap()
                    }
                    Array::Float(a) => {
                        PyArray1::from_vec(py, a.iter().filter_map(|r| r.ok()).collect())
                            .into_py_any(py)
                            .unwrap()
                    }
                    _ => py.None().into_py_any(py).unwrap(),
                },
                _ => py.None().into_py_any(py).unwrap(),
            };
            vec.push((key, obj));
        }
        vec
    }
}
