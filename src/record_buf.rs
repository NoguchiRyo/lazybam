use noodles::core::Position;
use noodles::sam::alignment::record::Flags;
use noodles::sam::alignment::record::MappingQuality;
use noodles::sam::alignment::record_buf::{Cigar, Data, QualityScores, Sequence};
use noodles::sam::alignment::RecordBuf;
use pyo3::prelude::*;

use crate::record_override;
#[pyclass]
pub struct PyRecordBuf {
    record_buf: RecordBuf,
}

impl PyRecordBuf {
    pub fn as_record_buf(&self) -> &RecordBuf {
        &self.record_buf
    }
}

#[pymethods]
impl PyRecordBuf {
    #[new]
    #[pyo3(signature = (qname, seq, qual, reference_sequence_id=None, cigar=None, alignment_start=None, mapping_quality=None, tags=None))]
    pub fn new(
        qname: String,
        seq: String,
        qual: Vec<u8>,
        reference_sequence_id: Option<u32>,
        cigar: Option<Vec<(u32, u32)>>,
        alignment_start: Option<u32>,
        mapping_quality: Option<u8>,
        tags: Option<Vec<(String, Py<PyAny>)>>,
    ) -> PyResult<Self> {
        let mut builder = RecordBuf::builder()
            .set_name(qname)
            .set_sequence(Sequence::from(seq.as_bytes()))
            .set_quality_scores(QualityScores::from(qual.clone()));

        let mut flag = Flags::UNMAPPED;

        if let Some(id) = reference_sequence_id {
            builder = builder.set_reference_sequence_id(id as usize);
            flag.remove(Flags::UNMAPPED);
        }
        if let Some(cigar) = cigar {
            let cigar: Cigar = record_override::convert_vec_to_cigar(cigar)
                .map_err(|e| PyErr::new::<pyo3::exceptions::PyValueError, _>(format!("{}", e)))?;
            builder = builder.set_cigar(cigar);
        }
        if let Some(start) = alignment_start {
            builder =
                builder.set_alignment_start(Position::try_from(start as usize).map_err(|e| {
                    PyErr::new::<pyo3::exceptions::PyValueError, _>(format!("{}", e))
                })?);
        }

        let mut tag_vec = Vec::new();
        if let Some(tag_list) = tags {
            for (k, v_any) in tag_list {
                if let (Ok(tag), Ok(val)) = (
                    record_override::convert_string_to_tag(k),
                    record_override::convert_pyany_to_value(v_any),
                ) {
                    tag_vec.push((tag, val));
                }
            }
        }
        if let Some((_tag, _value)) = tag_vec.first() {
            builder = builder.set_data(Data::from_iter(tag_vec.iter().cloned()));
        }

        if let Some(mq) = mapping_quality.and_then(MappingQuality::new) {
            builder = builder.set_mapping_quality(mq);
        }

        builder = builder.set_flags(flag);

        Ok(Self {
            record_buf: builder.build(),
        })
    }
}
