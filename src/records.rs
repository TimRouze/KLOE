use std::sync::Arc;

#[derive(Debug)]
pub struct SimplitigRecord {
    pub color_ids: Arc<Vec<u32>>,
    pub abundance_codes: Option<Arc<Vec<u8>>>,
    pub seq: Vec<u8>,
}

pub type SimplitigBatch = Vec<SimplitigRecord>;
