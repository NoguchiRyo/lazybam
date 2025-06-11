use noodles::sam::alignment::record::MappingQuality;
use noodles::sam::alignment::record_buf::Cigar;
use noodles::sam::alignment::record_buf::{QualityScores, Sequence as SeqBuf};
use noodles::sam::alignment::{
    record::cigar::op::Kind, record::cigar::Op, record::data::field::Tag,
    record_buf::data::field::Value,
};

use anyhow::Context;
use numpy::PyArrayMethods;
use numpy::{PyArray1, PyReadonlyArray1};
use pyo3::prelude::*;
use pyo3::types::PyAny;

/// Python 用に限定した「オーバーライド」構造体
#[pyclass]
#[derive(Clone)]
pub struct RecordOverride {
    pub qname: Option<String>,
    pub seq: Option<SeqBuf>,
    pub qual: Option<QualityScores>,
    pub reference_sequence_id: Option<u32>,
    pub cigar: Option<Cigar>,
    pub alignment_start: Option<u32>,
    pub tags: Vec<(Tag, Value)>,
    pub mapping_quality: Option<MappingQuality>,
}

#[pymethods]
impl RecordOverride {
    #[new]
    #[pyo3(signature = (qname=None, seq=None, qual=None, reference_sequence_id=None, cigar=None, alignment_start=None, tags=None, mapping_quality=None))]
    fn new(
        qname: Option<String>,
        seq: Option<String>,
        qual: Option<Vec<u8>>,
        reference_sequence_id: Option<u32>,
        cigar: Option<Vec<(u32, u32)>>,
        alignment_start: Option<u32>,
        tags: Option<Vec<(String, Py<PyAny>)>>,
        mapping_quality: Option<u8>,
    ) -> Self {
        let seq_opt = match seq {
            Some(s) => Some(SeqBuf::from(s.as_bytes())),
            None => None,
        };
        let qual_opt = match qual {
            Some(q) => Some(QualityScores::from(q)),
            None => None,
        };
        let cigar_opt = match cigar {
            Some(cigar_list) => convert_vec_to_cigar(cigar_list).ok(),
            None => None,
        };

        let mut tag_vec = Vec::new();
        if let Some(tag_list) = tags {
            for (k, v_any) in tag_list {
                if let (Ok(tag), Ok(val)) =
                    (convert_string_to_tag(k), convert_pyany_to_value(v_any))
                {
                    tag_vec.push((tag, val));
                }
            }
        }

        let mapq = match mapping_quality {
            Some(mq) => MappingQuality::new(mq),
            None => None,
        };

        RecordOverride {
            qname: qname,
            seq: seq_opt,
            qual: qual_opt,
            reference_sequence_id: reference_sequence_id,
            cigar: cigar_opt,
            alignment_start: alignment_start,
            tags: tag_vec,
            mapping_quality: mapq,
        }
    }

    /// override する reference_sequence_id (None なら元値を使う)
    #[setter]
    fn reference_sequence_id(&mut self, rid: u32) {
        self.reference_sequence_id = Some(rid);
    }

    #[setter]
    fn alignment_start(&mut self, pos: u32) {
        self.alignment_start = Some(pos);
    }

    #[setter]
    fn cigar(&mut self, cigar_list: Vec<(u32, u32)>) {
        // CIGAR の変換
        let cigar = convert_vec_to_cigar(cigar_list).unwrap();
        self.cigar = Some(cigar);
    }

    /// 追加タグ: Python からは List[(str, Any)] を受け取る
    #[setter]
    fn tags(&mut self, vals: Vec<(String, Py<PyAny>)>) {
        for (k, v_any) in vals {
            let tag = convert_string_to_tag(k).expect("Invalid tag");
            let val = convert_pyany_to_value(v_any).expect("Invalid value");
            self.tags.push((tag, val));
        }
    }
}

pub fn convert_string_to_tag(tag_str: String) -> anyhow::Result<Tag> {
    if tag_str.len() != 2 {
        return Err(anyhow::anyhow!("Invalid tag length: {}", tag_str.len()));
    }
    let tag_bytes = tag_str.as_bytes();
    let tag = Tag::try_from([tag_bytes[0], tag_bytes[1]])
        .map_err(|_| anyhow::anyhow!("Invalid tag bytes: {} {}", tag_bytes[0], tag_bytes[1]))?;
    Ok(tag)
}

pub fn convert_pyany_to_value(obj: PyObject) -> anyhow::Result<Value> {
    Python::with_gil(|py| {
        let any = obj.into_bound(py);
        if let Ok(i) = any.extract::<i64>() {
            // ← PyAnyMethods::extract :contentReference[oaicite:2]{index=2}
            return Ok(Value::try_from(i).context(format!("failed to convert Python int `{}`", i))?);
        }

        if let Ok(i) = any.extract::<i32>() {
            // ← PyAnyMethods::extract :contentReference[oaicite:2]{index=2}
            return Ok(Value::try_from(i).context(format!("failed to convert Python int `{}`", i))?);
        }

        if let Ok(i) = any.extract::<i16>() {
            // ← PyAnyMethods::extract :contentReference[oaicite:2]{index=2}
            return Ok(Value::try_from(i).context(format!("failed to convert Python int `{}`", i))?);
        }

        if let Ok(i) = any.extract::<i8>() {
            // ← PyAnyMethods::extract :contentReference[oaicite:2]{index=2}
            return Ok(Value::try_from(i).context(format!("failed to convert Python int `{}`", i))?);
        }

        if let Ok(i) = any.extract::<u32>() {
            // ← PyAnyMethods::extract :contentReference[oaicite:2]{index=2}
            return Ok(Value::try_from(i).context(format!("failed to convert Python int `{}`", i))?);
        }

        if let Ok(i) = any.extract::<u16>() {
            // ← PyAnyMethods::extract :contentReference[oaicite:2]{index=2}
            return Ok(Value::try_from(i).context(format!("failed to convert Python int `{}`", i))?);
        }

        if let Ok(i) = any.extract::<u8>() {
            // ← PyAnyMethods::extract :contentReference[oaicite:2]{index=2}
            return Ok(Value::try_from(i).context(format!("failed to convert Python int `{}`", i))?);
        }

        if let Ok(f) = any.extract::<f64>() {
            return Ok(Value::from(f as f32));
        }

        if let Ok(s) = any.extract::<String>() {
            return Ok(Value::from(s.as_str()));
        }

        // 2. If it is a 1D numpy array of i8, downcast and convert into Vec<i8>, then into BufValue::Array.
        if let Ok(py_arr) = any.downcast::<PyArray1<i8>>() {
            let array: PyReadonlyArray1<'_, i8> = py_arr.readonly();
            // as_slice() で &[i8] を取り出す
            let slice = array
                .as_slice()
                .map_err(|e| anyhow::anyhow!("Failed to get i8 slice: {}", e))?;
            let vec_i8 = slice.to_vec();
            // From<Vec<i8>> for BufValue があるので、これで配列タグが作れる
            let buf_value: Value = Vec::<i8>::from(vec_i8).into();
            return Ok(buf_value);
        }

        // 3. If it is a 1D numpy array of u8, downcast and convert into Vec<u8>, then into BufValue::Array.
        if let Ok(py_arr) = any.downcast::<PyArray1<u8>>() {
            let array: PyReadonlyArray1<'_, u8> = py_arr.readonly();
            let slice = array
                .as_slice()
                .map_err(|e| anyhow::anyhow!("Failed to get u8 slice: {}", e))?;
            let vec_u8 = slice.to_vec();
            let buf_value: Value = Vec::<u8>::from(vec_u8).into();
            return Ok(buf_value);
        }

        // 4. If it is a 1D numpy array of i16, downcast and convert into Vec<i16>, then into BufValue::Array.
        if let Ok(py_arr) = any.downcast::<PyArray1<i16>>() {
            let array: PyReadonlyArray1<'_, i16> = py_arr.readonly();
            let slice = array
                .as_slice()
                .map_err(|e| anyhow::anyhow!("Failed to get i16 slice: {}", e))?;
            let vec_i16 = slice.to_vec();
            let buf_value: Value = Vec::<i16>::from(vec_i16).into();
            return Ok(buf_value);
        }

        // 5. If it is a 1D numpy array of u16, downcast and convert into Vec<u16>, then into BufValue::Array.
        if let Ok(py_arr) = any.downcast::<PyArray1<u16>>() {
            let array: PyReadonlyArray1<'_, u16> = py_arr.readonly();
            let slice = array
                .as_slice()
                .map_err(|e| anyhow::anyhow!("Failed to get u16 slice: {}", e))?;
            let vec_u16 = slice.to_vec();
            let buf_value: Value = Vec::<u16>::from(vec_u16).into();
            return Ok(buf_value);
        }

        // 6. If it is a 1D numpy array of i32, downcast and convert into Vec<i32>, then into BufValue::Array.
        if let Ok(py_arr) = any.downcast::<PyArray1<i32>>() {
            let array: PyReadonlyArray1<'_, i32> = py_arr.readonly();
            let slice = array
                .as_slice()
                .map_err(|e| anyhow::anyhow!("Failed to get i32 slice: {}", e))?;
            let vec_i32 = slice.to_vec();
            let buf_value: Value = Vec::<i32>::from(vec_i32).into();
            return Ok(buf_value);
        }

        // 7. If it is a 1D numpy array of u32, downcast and convert into Vec<u32>, then into BufValue::Array.
        if let Ok(py_arr) = any.downcast::<PyArray1<u32>>() {
            let array: PyReadonlyArray1<'_, u32> = py_arr.readonly();
            let slice = array
                .as_slice()
                .map_err(|e| anyhow::anyhow!("Failed to get u32 slice: {}", e))?;
            let vec_u32 = slice.to_vec();
            let buf_value: Value = Vec::<u32>::from(vec_u32).into();
            return Ok(buf_value);
        }

        // 8. If it is a 1D numpy array of f32, downcast and convert into Vec<f32>, then into BufValue::Array.
        if let Ok(py_arr) = any.downcast::<PyArray1<f32>>() {
            let array: PyReadonlyArray1<'_, f32> = py_arr.readonly();
            let slice = array
                .as_slice()
                .map_err(|e| anyhow::anyhow!("Failed to get f32 slice: {}", e))?;
            let vec_f32 = slice.to_vec();
            let buf_value: Value = Vec::<f32>::from(vec_f32).into();
            return Ok(buf_value);
        }
        // その他はエラー
        let ty = any.get_type().name()?; // ← PyAnyMethods::get_type :contentReference[oaicite:3]{index=3}
        anyhow::bail!("unsupported Python type for Value: {}", ty);
    })
}

pub fn convert_vec_to_cigar(cigar_list: Vec<(u32, u32)>) -> anyhow::Result<Cigar> {
    let ops: Vec<Op> = cigar_list
        .into_iter()
        .map(|(k, l)| {
            let kind = match k {
                0 => Kind::Match,
                1 => Kind::SequenceMismatch,
                2 => Kind::Insertion,
                3 => Kind::Deletion,
                4 => Kind::Skip,
                5 => Kind::SoftClip,
                6 => Kind::HardClip,
                _ => return Err(anyhow::anyhow!("Invalid CIGAR operation: {}", k)),
            };
            Ok(Op::new(kind, l as usize))
        })
        .collect::<Result<_, _>>()?;
    Ok(Cigar::from(ops))
}
