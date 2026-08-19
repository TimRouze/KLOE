use core::panic;
use std::collections::HashMap;
use std::fs::{self, File};
use std::io::{self, BufRead, BufReader, BufWriter, Read, Result, Seek, Write};
use std::path::{Path, PathBuf};

use ggcat_api::{ExtraElaboration, GGCATConfig, GGCATInstance, GeneralSequenceBlockData};
use zstd::Decoder;

use crate::utils::vec2str;

const BUCKET_SIZES_MAGIC: &[u8; 4] = b"KSB2";
const POSITIONS_MAGIC: &[u8; 4] = b"KPS2";
const ID_TO_CID_MAGIC: &[u8; 4] = b"KIC2";

//   =========================================================================================== DECOMPRESSION ==============================================================================

/// High-level decompression entry point.
pub fn decompress(
    size_filename: &String,
    color_id_filename: &String,
    tigs_filename: &String,
    positions_filename: &String,
    filename_id: &String,
    out_dir: &String,
    wanted_files_path: &String,
    input_dir: String,
) -> std::io::Result<()> {
    decompress_with_options(
        size_filename,
        color_id_filename,
        tigs_filename,
        positions_filename,
        filename_id,
        out_dir,
        wanted_files_path,
        input_dir,
        GgcatRebuildConfig::default(),
    )
}

#[derive(Clone, Debug)]
pub struct GgcatRebuildConfig {
    pub enabled: bool,
    pub threads: usize,
    pub memory_gb: usize,
    pub k: usize,
    pub temp_dir: String,
    pub use_unitigs: bool,
    pub use_matchtigs: bool,
    pub use_eulertigs: bool,
}

impl Default for GgcatRebuildConfig {
    fn default() -> Self {
        Self {
            enabled: false,
            threads: 1,
            memory_gb: 8,
            k: 31,
            temp_dir: String::new(),
            use_unitigs: false,
            use_matchtigs: false,
            use_eulertigs: false,
        }
    }
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum RebuildMode {
    Simplitig,
    Unitig,
    Matchtig,
    Eulertig,
}

fn resolve_rebuild_mode(
    use_unitigs: bool,
    use_matchtigs: bool,
    use_eulertigs: bool,
) -> RebuildMode {
    if use_unitigs {
        RebuildMode::Unitig
    } else if use_matchtigs {
        RebuildMode::Matchtig
    } else if use_eulertigs {
        RebuildMode::Eulertig
    } else {
        RebuildMode::Simplitig
    }
}

fn ggcat_extra_elaboration(mode: RebuildMode) -> ExtraElaboration {
    match mode {
        RebuildMode::Unitig => ExtraElaboration::None,
        RebuildMode::Matchtig => ExtraElaboration::GreedyMatchtigs,
        RebuildMode::Eulertig => ExtraElaboration::FastEulertigs,
        RebuildMode::Simplitig => ExtraElaboration::FastSimplitigs,
    }
}

fn rebuild_mode_name(mode: RebuildMode) -> &'static str {
    match mode {
        RebuildMode::Simplitig => "simplitigs",
        RebuildMode::Unitig => "unitigs",
        RebuildMode::Matchtig => "matchtigs",
        RebuildMode::Eulertig => "eulertigs",
    }
}

fn normalize_out_dir(out_dir: &str) -> PathBuf {
    if out_dir.is_empty() {
        PathBuf::from(".")
    } else {
        PathBuf::from(out_dir)
    }
}

fn ensure_output_dir(out_dir: &str) -> io::Result<PathBuf> {
    let path = normalize_out_dir(out_dir);
    fs::create_dir_all(&path)?;
    Ok(path)
}

fn dump_output_path(out_dir: &str, source_path: &str) -> PathBuf {
    let out_dir_path = normalize_out_dir(out_dir);
    let trunc_filename = Path::new(source_path).file_stem().unwrap_or_default();
    out_dir_path.join(format!(
        "Dump_{}.fa",
        trunc_filename.to_str().unwrap_or("unknown")
    ))
}

fn collect_dump_fastas(out_dir: &str) -> std::io::Result<Vec<PathBuf>> {
    let out_dir_path = normalize_out_dir(out_dir);
    let mut dump_files = Vec::new();
    for entry_result in fs::read_dir(&out_dir_path)? {
        let entry = entry_result?;
        let path = entry.path();
        let is_dump = path
            .file_name()
            .and_then(|name| name.to_str())
            .map(|name| name.starts_with("Dump_") && name.ends_with(".fa"))
            .unwrap_or(false);
        if is_dump {
            dump_files.push(path);
        }
    }
    dump_files.sort();
    Ok(dump_files)
}

fn run_ggcat_rebuild(out_dir: &str, cfg: &GgcatRebuildConfig) -> std::io::Result<()> {
    let mode = resolve_rebuild_mode(cfg.use_unitigs, cfg.use_matchtigs, cfg.use_eulertigs);
    let dump_fastas = collect_dump_fastas(out_dir)?;
    if dump_fastas.is_empty() {
        return Err(io::Error::new(
            io::ErrorKind::NotFound,
            "no Dump_*.fa files found to rebuild with ggcat",
        ));
    }

    let out_dir_path = ensure_output_dir(out_dir)?;

    let rebuilt_output = out_dir_path.join(format!("rebuilt_{}.fa", rebuild_mode_name(mode)));
    let temp_dir = if cfg.temp_dir.is_empty() {
        out_dir_path.join("ggcat_rebuild_tmp")
    } else {
        PathBuf::from(&cfg.temp_dir)
    };
    fs::create_dir_all(&temp_dir)?;

    let ggcat_memory_gb = cfg.memory_gb.max(1);
    let instance = GGCATInstance::create(GGCATConfig {
        temp_dir: Some(temp_dir.clone()),
        memory: ggcat_memory_gb as f64,
        prefer_memory: false,
        total_threads_count: cfg.threads.max(1),
        intermediate_compression_level: None,
        stats_file: None,
        messages_callback: None,
    })
    .map_err(|err| io::Error::other(format!("create ggcat instance: {err}")))?;

    let streams = dump_fastas
        .iter()
        .map(|dump| {
            let resolved = fs::canonicalize(dump).unwrap_or_else(|_| dump.clone());
            GeneralSequenceBlockData::FASTA((resolved, None))
        })
        .collect::<Vec<_>>();

    println!(
        "Running embedded ggcat rebuild: output={}, mode={}, k={}, threads={}, memory={}GB",
        rebuilt_output.display(),
        rebuild_mode_name(mode),
        cfg.k,
        cfg.threads.max(1),
        ggcat_memory_gb
    );
    let graph_path = instance
        .build_graph(
            streams,
            rebuilt_output,
            None,
            None,
            cfg.k,
            cfg.threads.max(1),
            false,
            None,
            false,
            1,
            ggcat_extra_elaboration(mode),
            None,
        )
        .map_err(|err| io::Error::other(format!("run ggcat rebuild: {err}")))?;

    println!(
        "ggcat rebuild complete: {} (mode={})",
        graph_path.display(),
        rebuild_mode_name(mode)
    );
    Ok(())
}

pub fn decompress_with_options(
    size_filename: &String,
    color_id_filename: &String,
    tigs_filename: &String,
    positions_filename: &String,
    filename_id: &String,
    out_dir: &String,
    wanted_files_path: &String,
    input_dir: String,
    ggcat_cfg: GgcatRebuildConfig,
) -> std::io::Result<()> {
    println!("Writing decompressed data in {out_dir}");
    ensure_output_dir(out_dir)?;

    if wanted_files_path != "" {
        let input_file = File::open(input_dir.clone() + filename_id).unwrap();
        let input_reader = BufReader::new(input_file);
        let mut filenames_id_map = HashMap::new();
        let mut file_id: u32 = 0;
        for line_result in input_reader.lines() {
            let line = line_result?;
            if let Some((path, size)) = line.split_once(':') {
                filenames_id_map.insert(path.to_string(), (file_id, size.parse::<u64>().unwrap()));
                println!("File: {}, ID: {}, offset: {}", path, file_id, size);
            }
            file_id += 1;
        }

        let cid_map_out_filenames = match get_cid_to_id_targeted(
            &(input_dir.clone() + &color_id_filename),
            &filenames_id_map,
            &wanted_files_path,
        ) {
            Ok(map) => map,
            Err(e) => panic!("error gathering color set ids: {e:?}"),
        };
        let cid_to_id_map = cid_map_out_filenames.0;
        let wanted_filenames = cid_map_out_filenames.1;

        println!("Query file given, decompressing only subpart of archive....");
        decompress_wanted(
            &wanted_filenames,
            &(input_dir.clone() + positions_filename),
            cid_to_id_map,
            &(input_dir.to_owned() + tigs_filename),
            &(input_dir.to_owned() + size_filename),
            out_dir,
        );
    } else {
        println!("No query file given, decompressing entire archive....");
        let cid_to_id_map = match get_cid_to_id(&(input_dir.clone() + &color_id_filename)) {
            Ok(map) => map,
            Err(e) => panic!("Error getting cid to id map {e:?}"),
        };
        let input_file = File::open(input_dir.clone() + filename_id).unwrap();
        let input_reader = BufReader::new(input_file);
        let mut filenames_id = Vec::new();
        let mut file_id: u32 = 0;
        for line_result in input_reader.lines() {
            let line = line_result?;
            if let Some((path, _)) = line.split_once(":") {
                filenames_id.push((path.to_owned(), file_id));
            }
            file_id += 1;
        }
        println!("{}", input_dir.to_owned() + size_filename);
        decompress_all(
            &(input_dir.to_owned() + size_filename),
            &(input_dir.clone() + positions_filename),
            &(input_dir.to_owned() + tigs_filename),
            out_dir,
            filenames_id,
            cid_to_id_map,
        );
    }
    if ggcat_cfg.enabled {
        run_ggcat_rebuild(out_dir, &ggcat_cfg)?;
    }
    Ok(())
}

/// Preload all positions from the raw binary positions file.
///
/// The positions file contains N entries of exactly 16 bytes each:
/// [u64 tigs_pos LE][u64 sizes_pos LE].
/// Returns a Vec where index i corresponds to CID i.
fn read_varint_u64_from_reader(reader: &mut impl Read) -> io::Result<Option<u64>> {
    let mut shift = 0u32;
    let mut value = 0u64;
    let mut saw_byte = false;
    loop {
        let mut b = [0u8; 1];
        match reader.read_exact(&mut b) {
            Ok(()) => {}
            Err(err) if err.kind() == io::ErrorKind::UnexpectedEof && !saw_byte => return Ok(None),
            Err(err) if err.kind() == io::ErrorKind::UnexpectedEof => {
                return Err(io::Error::new(
                    io::ErrorKind::UnexpectedEof,
                    "truncated varint value",
                ))
            }
            Err(err) => return Err(err),
        }
        saw_byte = true;
        let byte = b[0];
        value |= ((byte & 0x7f) as u64) << shift;
        if byte & 0x80 == 0 {
            return Ok(Some(value));
        }
        shift += 7;
        if shift >= 64 {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "varint overflow while decoding u64",
            ));
        }
    }
}

fn decode_varint_delta_cids(payload: &[u8], context: &str) -> Result<Vec<usize>> {
    let mut cursor = std::io::Cursor::new(payload);
    let mut cids = Vec::new();
    let mut current = 0usize;
    loop {
        let Some(delta_u64) = read_varint_u64_from_reader(&mut cursor)? else {
            break;
        };
        let delta = usize::try_from(delta_u64).map_err(|_| {
            io::Error::new(
                io::ErrorKind::InvalidData,
                format!(
                    "cid delta {} cannot fit usize while decoding '{}'",
                    delta_u64, context
                ),
            )
        })?;
        current = current.checked_add(delta).ok_or_else(|| {
            io::Error::new(
                io::ErrorKind::InvalidData,
                format!("cid delta overflow while decoding '{}'", context),
            )
        })?;
        cids.push(current);
    }
    Ok(cids)
}

fn decode_varint_deltas_to_sizes(payload: &[u8], context: &str) -> Result<Vec<usize>> {
    let mut cursor = std::io::Cursor::new(payload);
    let mut sizes = Vec::new();
    let mut prev = 0usize;
    loop {
        let Some(delta_u64) = read_varint_u64_from_reader(&mut cursor)? else {
            break;
        };
        let delta = usize::try_from(delta_u64).map_err(|_| {
            io::Error::new(
                io::ErrorKind::InvalidData,
                format!(
                    "size delta {} cannot fit usize while decoding '{}'",
                    delta_u64, context
                ),
            )
        })?;
        let size = prev.checked_add(delta).ok_or_else(|| {
            io::Error::new(
                io::ErrorKind::InvalidData,
                format!("size delta overflow while decoding '{}'", context),
            )
        })?;
        sizes.push(size);
        prev = size;
    }
    Ok(sizes)
}

fn preload_positions(positions_filename: &str) -> Result<Vec<(u64, u64)>> {
    let mut file = BufReader::new(File::open(positions_filename)?);
    let mut magic = [0u8; 4];
    match file.read_exact(&mut magic) {
        Ok(()) => {}
        Err(err) if err.kind() == io::ErrorKind::UnexpectedEof => return Ok(Vec::new()),
        Err(err) => return Err(err),
    }

    if &magic == POSITIONS_MAGIC {
        let entries = read_varint_u64_from_reader(&mut file)?.ok_or_else(|| {
            io::Error::new(
                io::ErrorKind::UnexpectedEof,
                "missing positions entry count in compact positions file",
            )
        })?;
        let entries = usize::try_from(entries).map_err(|_| {
            io::Error::new(
                io::ErrorKind::InvalidData,
                "positions entry count cannot fit usize",
            )
        })?;
        let mut positions = Vec::with_capacity(entries);
        let mut tigs_pos = 0u64;
        let mut sizes_pos = 0u64;
        for _ in 0..entries {
            let dt = read_varint_u64_from_reader(&mut file)?.ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::UnexpectedEof,
                    "truncated tigs delta in positions file",
                )
            })?;
            let ds = read_varint_u64_from_reader(&mut file)?.ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::UnexpectedEof,
                    "truncated sizes delta in positions file",
                )
            })?;
            tigs_pos = tigs_pos.checked_add(dt).ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    "tigs position delta overflow in positions file",
                )
            })?;
            sizes_pos = sizes_pos.checked_add(ds).ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    "sizes position delta overflow in positions file",
                )
            })?;
            positions.push((tigs_pos, sizes_pos));
        }
        println!(
            "Preloaded {} compact position entries from {}",
            positions.len(),
            positions_filename
        );
        return Ok(positions);
    }

    file.seek(std::io::SeekFrom::Start(0))?;
    let file_size = file.get_ref().metadata()?.len() as usize;
    let num_entries = file_size / 16;
    let mut positions = Vec::with_capacity(num_entries);
    let mut buf = [0u8; 16];
    for _ in 0..num_entries {
        file.read_exact(&mut buf)?;
        let tigs_pos = u64::from_le_bytes(buf[..8].try_into().unwrap());
        let sizes_pos = u64::from_le_bytes(buf[8..16].try_into().unwrap());
        positions.push((tigs_pos, sizes_pos));
    }
    println!(
        "Preloaded {} legacy position entries from {}",
        num_entries, positions_filename
    );
    Ok(positions)
}

/// Decode comma-separated delta-encoded CIDs into absolute CID indexes.
///
/// Example: "3,2,0,5" => [3,5,5,10]
fn decode_delta_cids(text: &str, context: &str) -> Result<Vec<usize>> {
    let mut cids = Vec::new();
    let mut current = 0usize;
    for token in text.split(',') {
        let token = token.trim();
        if token.is_empty() {
            continue;
        }
        let delta = token.parse::<usize>().map_err(|err| {
            io::Error::new(
                io::ErrorKind::InvalidData,
                format!("invalid cid delta '{}' in '{}': {}", token, context, err),
            )
        })?;
        current = current.checked_add(delta).ok_or_else(|| {
            io::Error::new(
                io::ErrorKind::InvalidData,
                format!("cid delta overflow while decoding '{}'", context),
            )
        })?;
        cids.push(current);
    }
    Ok(cids)
}

fn decode_cid_entry_payload(
    payload: &[u8],
    context: &str,
    binary_varints: bool,
) -> Result<Vec<usize>> {
    if binary_varints {
        decode_varint_delta_cids(payload, context)
    } else {
        let text = String::from_utf8(payload.to_vec()).map_err(|err| {
            io::Error::new(
                io::ErrorKind::InvalidData,
                format!("invalid UTF-8 while decoding '{}': {}", context, err),
            )
        })?;
        decode_delta_cids(&text, context)
    }
}

/// Preload all bucket sizes from the sizes file into a Vec indexed by byte offset.
///
/// The sizes file format: repeated [u64 compressed_len][compressed_data...].
/// Returns a HashMap from byte_offset -> Vec<usize> of actual sizes.
fn preload_sizes(size_filename: &str) -> Result<HashMap<u64, Vec<usize>>> {
    let mut file = BufReader::new(File::open(size_filename)?);
    let mut sizes_map = HashMap::new();
    let mut magic = [0u8; 4];
    match file.read_exact(&mut magic) {
        Ok(()) => {}
        Err(err) if err.kind() == io::ErrorKind::UnexpectedEof => return Ok(HashMap::new()),
        Err(err) => return Err(err),
    }

    if &magic == BUCKET_SIZES_MAGIC {
        let mut global_group_index = 0u64;
        loop {
            let mut groups_buf = [0u8; 4];
            match file.read_exact(&mut groups_buf) {
                Ok(()) => {}
                Err(err) if err.kind() == io::ErrorKind::UnexpectedEof => break,
                Err(err) => return Err(err),
            }
            let groups = u32::from_le_bytes(groups_buf) as usize;
            if groups == 0 {
                break;
            }
            let mut compressed_len_buf = [0u8; 8];
            file.read_exact(&mut compressed_len_buf)?;
            let compressed_len = u64::from_le_bytes(compressed_len_buf) as usize;

            let offsets_len = groups + 1;
            let mut offsets = vec![0u32; offsets_len];
            for off in &mut offsets {
                let mut buf = [0u8; 4];
                file.read_exact(&mut buf)?;
                *off = u32::from_le_bytes(buf);
            }

            let mut compressed = vec![0u8; compressed_len];
            file.read_exact(&mut compressed)?;

            let mut decompressed = Vec::new();
            Decoder::new(&compressed[..])?.read_to_end(&mut decompressed)?;

            for i in 0..groups {
                let start = offsets[i] as usize;
                let end = offsets[i + 1] as usize;
                if end < start || end > decompressed.len() {
                    return Err(io::Error::new(
                        io::ErrorKind::InvalidData,
                        "invalid bucket block offsets while preloading sizes",
                    ));
                }
                let key = global_group_index + i as u64;
                let context = format!("{}#group{}", size_filename, key);
                let sizes = decode_varint_deltas_to_sizes(&decompressed[start..end], &context)?;
                sizes_map.insert(key, sizes);
            }
            global_group_index += groups as u64;
        }

        println!(
            "Preloaded {} compact size buckets from {}",
            sizes_map.len(),
            size_filename
        );
        return Ok(sizes_map);
    }

    file.seek(std::io::SeekFrom::Start(0))?;
    let mut offset: u64 = 0;
    loop {
        let mut len_buf = [0u8; 8];
        match file.read_exact(&mut len_buf) {
            Ok(()) => {}
            Err(e) if e.kind() == io::ErrorKind::UnexpectedEof => break,
            Err(e) => return Err(e),
        }
        let compressed_size = u64::from_le_bytes(len_buf) as usize;
        if compressed_size == 0 {
            break;
        }
        let mut compressed_buffer = vec![0; compressed_size];
        file.read_exact(&mut compressed_buffer)?;

        let mut decoder = Decoder::new(&compressed_buffer[..])?;
        let mut decompressed = Vec::new();
        decoder.read_to_end(&mut decompressed)?;

        let mut sizes = Vec::new();
        let mut prev = 0usize;
        for chunk in decompressed.chunks_exact(8) {
            let delta = usize::from_le_bytes(chunk.try_into().unwrap());
            let actual_size = prev
                .checked_add(delta)
                .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidData, "size overflow"))?;
            sizes.push(actual_size);
            prev = actual_size;
        }

        sizes_map.insert(offset, sizes);
        offset += 8 + compressed_size as u64;
    }
    println!(
        "Preloaded {} legacy size buckets from {}",
        sizes_map.len(),
        size_filename
    );
    Ok(sizes_map)
}

/// Decompress the entire archive to individual FASTA shards.
///
/// Uses preloaded positions and sizes to avoid per-CID file reopens.
/// Keeps output file handles open in a HashMap to avoid per-unitig open/close.
fn decompress_all(
    size_filename: &String,
    positions_filename: &String,
    tigs_filename: &String,
    out_dir: &String,
    filenames: Vec<(String, u32)>,
    cid_to_id_map: HashMap<usize, Vec<u32>>,
) {
    let mut tigs_file =
        BufReader::new(File::open(&tigs_filename).expect("Error opening tigs file"));

    let total_cids = cid_to_id_map.len();
    println!("nb cid: {}", total_cids);

    // Preload all positions and sizes into memory
    let all_positions = preload_positions(positions_filename).expect("Failed to preload positions");
    let all_sizes = preload_sizes(size_filename).expect("Failed to preload sizes");

    // Keep output file handles open
    let mut writers: HashMap<u32, BufWriter<File>> = HashMap::new();

    let mut sorted_cids: Vec<_> = cid_to_id_map.keys().cloned().collect();
    sorted_cids.sort();

    let mut total_unitigs = 0_u64;
    for (idx, cid) in sorted_cids.iter().enumerate() {
        let file_ids = cid_to_id_map.get(cid).unwrap();

        // Progress logging
        if idx % 1000 == 0 || idx == total_cids - 1 {
            println!(
                "Processing CID {}/{} ({} unitigs written so far)",
                idx + 1,
                total_cids,
                total_unitigs
            );
        }

        // Look up position from preloaded data
        let pos_index = *cid;
        if pos_index >= all_positions.len() {
            eprintln!(
                "Error: position index {} out of range for CID {}",
                pos_index, cid
            );
            continue;
        }
        let (tigs_pos, sizes_pos) = all_positions[pos_index];

        // Look up sizes from preloaded data
        let sizes = match all_sizes.get(&sizes_pos) {
            Some(s) => s,
            None => {
                eprintln!(
                    "Error: no sizes found at offset {} for CID {}",
                    sizes_pos, cid
                );
                continue;
            }
        };

        if sizes.is_empty() {
            continue;
        }

        tigs_file
            .seek(std::io::SeekFrom::Start(tigs_pos))
            .expect("Failed to seek in tigs file");

        for size in sizes {
            if *size < 31 {
                eprintln!("Warning: unitig size {} is less than k-mer size", size);
                continue;
            }

            let read_size = size.div_ceil(4);
            let mut tig_buffer = vec![0; read_size];
            tigs_file
                .read_exact(&mut tig_buffer)
                .expect("Failed to read tig");

            let tig = vec2str(&tig_buffer, size);
            for file_id in file_ids {
                let writer = writers.entry(*file_id).or_insert_with(|| {
                    let curr_filename = &filenames[*file_id as usize];
                    let output_path = dump_output_path(out_dir, &curr_filename.0);
                    BufWriter::new(
                        File::options()
                            .append(true)
                            .create(true)
                            .open(output_path)
                            .expect("Unable to create file"),
                    )
                });
                writeln!(writer, ">").unwrap();
                writeln!(writer, "{}", tig).unwrap();
            }
            total_unitigs += 1;
        }
    }

    // Flush all writers at the end
    for (_, mut writer) in writers {
        writer.flush().unwrap();
    }
    println!(
        "Decompression complete: {} unitigs written across {} CIDs",
        total_unitigs, total_cids
    );
}

/// Read full id->color_id file and build color id -> list of file ids.
fn get_cid_to_id(color_id_filename: &String) -> Result<HashMap<usize, Vec<u32>>> {
    let mut color_id_file = BufReader::new(
        File::open(color_id_filename)
            .expect("Error opening color id file, are you sure you gave the right path?"),
    );
    let mut cid_ids_map = HashMap::new();
    let mut counter: u32 = 0;
    println!("Reading CID to ID file: {}", color_id_filename);

    let mut magic = [0u8; 4];
    let binary_varints = match color_id_file.read_exact(&mut magic) {
        Ok(()) if &magic == ID_TO_CID_MAGIC => true,
        Ok(()) => {
            color_id_file.seek(std::io::SeekFrom::Start(0))?;
            false
        }
        Err(err) if err.kind() == io::ErrorKind::UnexpectedEof => return Ok(HashMap::new()),
        Err(err) => return Err(err),
    };

    let mut buffer_size = [0; 8];
    color_id_file.read_exact(&mut buffer_size)?;
    let mut size_read = u64::from_le_bytes(buffer_size) as usize;
    while size_read != 0 {
        let mut buffer = vec![0; size_read];
        color_id_file.read_exact(&mut buffer)?;
        let mut decompressed_data = Vec::new();
        {
            let mut decoder = Decoder::new(&buffer[..])?;
            decoder.read_to_end(&mut decompressed_data)?;
        }
        let context = format!("{}#{}", color_id_filename, counter);
        let decoded_cids = decode_cid_entry_payload(&decompressed_data, &context, binary_varints)?;
        for cid in decoded_cids {
            cid_ids_map
                .entry(cid)
                .and_modify(|list: &mut Vec<_>| list.push(counter))
                .or_insert(Vec::from([counter]));
        }
        color_id_file.read_exact(&mut buffer_size)?;
        size_read = u64::from_le_bytes(buffer_size) as usize;
        counter += 1;
    }

    Ok(cid_ids_map)
}

/// Build cid -> file id mapping for a targeted subset of files.
fn get_cid_to_id_targeted(
    color_id_filename: &String,
    filenames_id_map: &HashMap<String, (u32, u64)>,
    wanted_files_path: &String,
) -> std::io::Result<(HashMap<usize, Vec<u32>>, Vec<(String, u32)>)> {
    let mut color_id_file = BufReader::new(
        File::open(color_id_filename)
            .expect("Error opening color id file, are you sure you gave the right path?"),
    );
    let mut cid_ids_map = HashMap::new();
    let mut wanted_filenames = Vec::new();

    let mut magic = [0u8; 4];
    let binary_varints = match color_id_file.read_exact(&mut magic) {
        Ok(()) if &magic == ID_TO_CID_MAGIC => true,
        Ok(()) => {
            color_id_file.seek(std::io::SeekFrom::Start(0))?;
            false
        }
        Err(err) if err.kind() == io::ErrorKind::UnexpectedEof => {
            return Ok((HashMap::new(), Vec::new()))
        }
        Err(err) => return Err(err),
    };

    let wanted_file = File::open(wanted_files_path)?;
    let wanted_reader = BufReader::new(wanted_file);
    for line_result in wanted_reader.lines() {
        let line = line_result?;
        if filenames_id_map.contains_key(&line) {
            let entry = filenames_id_map.get(&line).unwrap();
            color_id_file.seek(std::io::SeekFrom::Start(entry.1))?;
            let mut buffer_size = [0; 8];
            color_id_file.read_exact(&mut buffer_size)?;
            let size_read = u64::from_le_bytes(buffer_size) as usize;
            let mut buffer = vec![0; size_read];
            color_id_file.read_exact(&mut buffer)?;
            let mut decompressed_data = Vec::new();
            {
                let mut decoder = Decoder::new(&buffer[..])?;
                decoder.read_to_end(&mut decompressed_data)?;
            }
            let context = format!("{}@{}", color_id_filename, entry.1);
            let decoded_cids =
                decode_cid_entry_payload(&decompressed_data, &context, binary_varints)?;
            for cid in decoded_cids {
                cid_ids_map
                    .entry(cid)
                    .and_modify(|list: &mut Vec<_>| list.push(entry.0))
                    .or_insert(Vec::from([entry.0]));
            }
            wanted_filenames.push((line.clone(), entry.0));
        } else {
            println!(
                "FILE {} NOT FOUND IN ARCHIVE, CHECK SPELLING OR ACTUAL PRESENCE IN ARCHIVE",
                line
            );
        }
    }
    Ok((cid_ids_map, wanted_filenames))
}

/// Decompress only the wanted files from the archive.
///
/// Uses preloaded positions and sizes. Keeps output file handles open.
fn decompress_wanted(
    wanted_files: &Vec<(String, u32)>,
    positions_filename: &String,
    cid_to_id_map: HashMap<usize, Vec<u32>>,
    tigs_filename: &String,
    size_filename: &String,
    out_dir: &String,
) {
    let total_cids = cid_to_id_map.len();
    println!("NB COLOR TO DECOMPRESS: {}", total_cids);
    println!("Wanted files:");
    for elem in wanted_files {
        println!("{} : {}", elem.0, elem.1);
    }

    let mut tigs_file =
        BufReader::new(File::open(&tigs_filename).expect("Error opening tigs file"));

    // Preload all positions and sizes into memory
    let all_positions = preload_positions(positions_filename).expect("Failed to preload positions");
    let all_sizes = preload_sizes(size_filename).expect("Failed to preload sizes");

    // Build a set of wanted file_ids for quick lookup
    let wanted_ids: std::collections::HashSet<u32> = wanted_files.iter().map(|w| w.1).collect();

    // Keep output file handles open
    let mut writers: HashMap<u32, BufWriter<File>> = HashMap::new();
    // Pre-open all wanted output files
    for wanted_file in wanted_files {
        let output_path = dump_output_path(out_dir, &wanted_file.0);
        let writer = BufWriter::new(
            File::options()
                .append(true)
                .create(true)
                .open(output_path)
                .expect("Unable to create file"),
        );
        writers.insert(wanted_file.1, writer);
    }

    let mut sorted_cids: Vec<_> = cid_to_id_map.keys().cloned().collect();
    sorted_cids.sort();

    let mut total_unitigs = 0_u64;
    for (idx, cid) in sorted_cids.iter().enumerate() {
        let file_ids = cid_to_id_map.get(cid).unwrap();

        // Progress logging
        if idx % 1000 == 0 || idx == total_cids - 1 {
            println!(
                "Processing CID {}/{} ({} unitigs written so far)",
                idx + 1,
                total_cids,
                total_unitigs
            );
        }

        // Look up position from preloaded data
        let pos_index = *cid;
        if pos_index >= all_positions.len() {
            eprintln!(
                "Error: position index {} out of range for CID {}",
                pos_index, cid
            );
            continue;
        }
        let (tigs_pos, sizes_pos) = all_positions[pos_index];

        // Look up sizes from preloaded data
        let sizes = match all_sizes.get(&sizes_pos) {
            Some(s) => s,
            None => {
                eprintln!(
                    "Error: no sizes found at offset {} for CID {}",
                    sizes_pos, cid
                );
                continue;
            }
        };

        if sizes.is_empty() {
            continue;
        }
        tigs_file
            .seek(std::io::SeekFrom::Start(tigs_pos))
            .expect("Failed to seek in tigs file");

        for size in sizes {
            if *size < 31 {
                eprintln!("Warning: unitig size {} is less than k-mer size", size);
                continue;
            }

            let read_size = size.div_ceil(4);
            let mut tig_buffer = vec![0; read_size];
            tigs_file
                .read_exact(&mut tig_buffer)
                .expect("Failed to read tig");

            let tig = vec2str(&tig_buffer, size);

            for file_id in file_ids {
                if wanted_ids.contains(file_id) {
                    if let Some(writer) = writers.get_mut(file_id) {
                        writeln!(writer, ">").unwrap();
                        writeln!(writer, "{}", tig).unwrap();
                    }
                }
            }
            total_unitigs += 1;
        }
    }

    // Flush all writers at the end
    for (_, mut writer) in writers {
        writer.flush().unwrap();
    }
    println!(
        "Decompression complete: {} unitigs written across {} CIDs",
        total_unitigs, total_cids
    );
}
