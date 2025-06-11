use numpy::PyArray1;
use pyo3::prelude::*;
use pyo3::IntoPyObjectExt;

use noodles::sam::alignment::record::data::field::value::Array;
use noodles::sam::alignment::record::data::field::Tag;
use noodles::sam::alignment::record_buf::{Cigar, Data, QualityScores, Sequence as SeqBuf};
use noodles::sam::alignment::{
    record::Flags, //, MappingQuality
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
        let mut qname_opt = self.qname();
        let mut seq_opt = SeqBuf::from(self.seq().into_bytes());
        let mut qual_opt = QualityScores::from(self.record.quality_scores().as_ref().to_vec());
        let mut data = Data::try_from(self.record.data()).unwrap_or_default();
        let mut mapq_opt = self.record.mapping_quality();

        let mut position_opt = match self.record.alignment_start() {
            Some(Ok(pos)) => Some(pos),
            Some(Err(_)) => return Err(anyhow::anyhow!("Invalid alignment start position")),
            None => None,
        };

        let mut ref_id_opt = match self.record.reference_sequence_id() {
            Some(Ok(rid)) => Some(rid),
            Some(Err(_)) => return Err(anyhow::anyhow!("Invalid reference sequence ID")),
            None => None,
        };

        let mut cigar_vec: Vec<Op> = self.record.cigar().iter().filter_map(Result::ok).collect();

        let mut flag = self.record.flags();

        if let Some(ov) = &self.record_override {
            for (tag, value) in &ov.tags {
                data.insert(*tag, value.clone());
            }
            if let Some(cigar) = &ov.cigar {
                cigar_vec = cigar.iter().filter_map(Result::ok).collect();
            }
            if let Some(rid) = ov.reference_sequence_id {
                ref_id_opt = Some(rid as usize);
                flag.remove(Flags::UNMAPPED);
            }
            if let Some(start) = ov.alignment_start {
                position_opt = Some(Position::try_from(start as usize)?);
            }
            if let Some(seq) = &ov.seq {
                seq_opt = SeqBuf::from(seq.clone());
            }
            if let Some(qual) = &ov.qual {
                qual_opt = QualityScores::from(qual.clone());
            }
            if let Some(qname) = &ov.qname {
                qname_opt = qname.clone();
            }
            if let Some(mapq) = &ov.mapping_quality {
                mapq_opt = Some(*mapq)
            }
        }
        // builder
        let mut builder = RecordBuf::builder()
            .set_name(qname_opt)
            .set_flags(flag)
            .set_cigar(Cigar::from(cigar_vec))
            .set_sequence(seq_opt)
            .set_quality_scores(qual_opt)
            .set_data(data);

        if let Some(rid) = ref_id_opt {
            builder = builder.set_reference_sequence_id(rid);
        }
        if let Some(pos) = position_opt {
            builder = builder.set_alignment_start(pos);
        }
        if let Some(mapq) = mapq_opt {
            builder = builder.set_mapping_quality(mapq);
        }
        let record_buf = builder.build();

        Ok(record_buf)
    }
}

#[pymethods]
impl PyBamRecord {
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
    fn rid(&self) -> i32 {
        self.record
            .reference_sequence_id()
            .and_then(|r| r.ok())
            .map(|r| r as i32)
            .unwrap_or(-1)
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

    fn get_field_by_tag<'py>(&self, tag: &str, py: Python<'py>) -> PyResult<PyObject> {
        // First, convert tag to two bytes
        let tag_bytes = tag.as_bytes();
        // tag が 2 バイトでない場合はエラー
        if tag_bytes.len() != 2 {
            return Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(
                "tag must be 2 bytes",
            ));
        }
        // それ以外は元の record.data() から取得
        for result in self.record.data().iter() {
            let (key, value) = result.map_err(|_| {
                PyErr::new::<pyo3::exceptions::PyValueError, _>(format!(
                    "failed to get field by tag: {}",
                    tag
                ))
            })?;
            if key == Tag::new(tag_bytes[0], tag_bytes[1]) {
                let obj = match value {
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
                return Ok(obj);
            }
        }
        Err(PyErr::new::<pyo3::exceptions::PyKeyError, _>(format!(
            "tag not found: {}",
            tag
        )))
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
