use crate::memory_fs::file::internal::{FileChunk, MemoryFileInternal, OpenMode};
use parking_lot::lock_api::ArcRwLockReadGuard;
use parking_lot::{RawRwLock, RwLock};
use std::cmp::min;
use std::io;
use std::io::{ErrorKind, Read, Seek, SeekFrom};
use std::ops::Range;
use std::path::{Path, PathBuf};
use std::slice::from_raw_parts;
use std::sync::Arc;

use super::writer::FileWriter;

#[derive(Clone)]
pub struct FileRangeReference {
    file: Arc<RwLock<MemoryFileInternal>>,
    start_chunk: usize,
    start_chunk_offset: usize,
    bytes_count: usize,
}

impl FileRangeReference {
    pub unsafe fn copy_to_unsync(&self, other: &FileWriter) {
        let file = self.file.read();
        let mut chunk_index = self.start_chunk;
        let mut chunk_offset = self.start_chunk_offset;
        let mut written_bytes = 0;

        while written_bytes < self.bytes_count {
            let chunk = file.get_chunk(chunk_index);
            let chunk = chunk.read();
            let to_copy = (chunk.get_length() - chunk_offset).min(self.bytes_count - written_bytes);

            let data = chunk
                .get_ptr(&file.get_underlying_file(), None)
                .add(chunk_offset);

            other.write_all_parallel(from_raw_parts(data, to_copy), to_copy);

            written_bytes += to_copy;
            chunk_index += 1;
            chunk_offset = 0;
        }

        // other.write_all_parallel(data, el_size)
    }
}

pub struct FileReader {
    path: PathBuf,
    file: Arc<RwLock<MemoryFileInternal>>,
    current_chunk_ref: Option<ArcRwLockReadGuard<RawRwLock, FileChunk>>,
    current_chunk_index: usize,
    chunks_count: usize,
    current_ptr: *const u8,
    current_len: usize,
    current_file_position: usize,
    prefetch_amount: Option<usize>,
}

unsafe impl Sync for FileReader {}
unsafe impl Send for FileReader {}

impl Clone for FileReader {
    fn clone(&self) -> Self {
        Self {
            path: self.path.clone(),
            file: self.file.clone(),
            current_chunk_ref: self
                .current_chunk_ref
                .as_ref()
                .map(|c| ArcRwLockReadGuard::rwlock(&c).read_arc()),
            current_chunk_index: self.current_chunk_index,
            chunks_count: self.chunks_count,
            current_ptr: self.current_ptr,
            current_len: self.current_len,
            current_file_position: self.current_file_position,
            prefetch_amount: self.prefetch_amount,
        }
    }
}

impl FileReader {
    fn set_chunk_info(&mut self, index: usize) {
        let file = self.file.read();

        let chunk = file.get_chunk(index);
        let chunk_guard = chunk.read_arc();

        let underlying_file = file.get_underlying_file();

        self.current_ptr = chunk_guard.get_ptr(&underlying_file, self.prefetch_amount);
        self.current_len = chunk_guard.get_length();
        self.current_chunk_ref = Some(chunk_guard);
    }

    pub fn open(path: impl AsRef<Path>, prefetch_amount: Option<usize>) -> Option<Self> {
        let file = match MemoryFileInternal::retrieve_reference(&path) {
            None => MemoryFileInternal::create_from_fs(&path)?,
            Some(x) => x,
        };

        let mut file_lock = file.write();

        file_lock.open(OpenMode::Read).unwrap();
        let chunks_count = file_lock.get_chunks_count();
        drop(file_lock);

        let mut reader = Self {
            path: path.as_ref().into(),
            file,
            current_chunk_ref: None,
            current_chunk_index: 0,
            chunks_count,
            current_ptr: std::ptr::null(),
            current_len: 0,
            current_file_position: 0,
            prefetch_amount,
        };

        if reader.chunks_count > 0 {
            reader.set_chunk_info(0);
        }

        Some(reader)
    }

    pub fn get_unique_file_id(&self) -> usize {
        self.file.data_ptr() as usize
    }

    pub fn total_file_size(&self) -> usize {
        self.file.read().len()
    }

    pub fn get_file_path(&self) -> &Path {
        &self.path
    }

    pub fn close_and_remove(self, remove_fs: bool) -> bool {
        MemoryFileInternal::delete(self.path, remove_fs)
    }

    pub fn get_range_reference(&self, file_range: Range<u64>) -> FileRangeReference {
        let file = self.file.read();
        let mut chunk_index = 0;
        let mut chunk_offset = 0;
        let mut start = file_range.start;

        while start > 0 {
            let chunk = file.get_chunk(chunk_index);
            let chunk = chunk.read();
            let len = chunk.get_length() as u64;
            if start < len {
                chunk_offset = start as usize;
                break;
            }
            start -= len;
            chunk_index += 1;
        }

        FileRangeReference {
            file: self.file.clone(),
            start_chunk: chunk_index,
            start_chunk_offset: chunk_offset,
            bytes_count: file_range.end as usize - file_range.start as usize,
        }
    }

    // pub fn get_typed_chunks_mut<T>(&mut self) -> Option<impl Iterator<Item = &mut [T]>> {
    //     todo!();
    //     Some((0..1).into_iter().map(|_| &mut [0, 1][..]))
    // }
}

impl Read for FileReader {
    fn read(&mut self, buf: &mut [u8]) -> std::io::Result<usize> {
        let mut bytes_written = 0;

        while bytes_written != buf.len() {
            if self.current_len == 0 {
                self.current_chunk_index += 1;
                if self.current_chunk_index >= self.chunks_count {
                    // End of file
                    return Ok(bytes_written);
                }
                self.set_chunk_info(self.current_chunk_index);
            }

            let copyable_bytes = min(buf.len() - bytes_written, self.current_len);

            unsafe {
                std::ptr::copy_nonoverlapping(
                    self.current_ptr,
                    buf.as_mut_ptr().add(bytes_written),
                    copyable_bytes,
                );
                self.current_ptr = self.current_ptr.add(copyable_bytes);
                self.current_len -= copyable_bytes;
                self.current_file_position += copyable_bytes;
                bytes_written += copyable_bytes;
            }
        }

        Ok(bytes_written)
    }

    #[inline(always)]
    fn read_exact(&mut self, buf: &mut [u8]) -> std::io::Result<()> {
        match self.read(buf) {
            Ok(count) => {
                if count == buf.len() {
                    Ok(())
                } else {
                    Err(io::Error::new(
                        io::ErrorKind::Other,
                        "Unexpected error while reading",
                    ))
                }
            }
            Err(err) => Err(err),
        }
    }
}

impl Seek for FileReader {
    fn seek(&mut self, pos: SeekFrom) -> io::Result<u64> {
        match pos {
            SeekFrom::Start(mut offset) => {
                let mut chunk_idx = 0;

                let file = self.file.read();
                while chunk_idx < self.chunks_count {
                    let len = file.get_chunk(chunk_idx).read().get_length();
                    if offset < (len as u64) {
                        break;
                    }
                    chunk_idx += 1;
                    offset -= len as u64;
                }

                if chunk_idx == self.chunks_count {
                    return Err(std::io::Error::new(
                        ErrorKind::UnexpectedEof,
                        "Unexpected eof",
                    ));
                }

                self.current_chunk_index = chunk_idx;
                drop(file);
                self.set_chunk_info(chunk_idx);

                self.current_ptr = unsafe { self.current_ptr.add(offset as usize) };
                self.current_len -= offset as usize;
                self.current_file_position = offset as usize;

                return Ok(offset);
            }
            SeekFrom::Current(offset) => {
                assert!(offset >= 0); // < 0  not supported yet
                let mut offset = offset as usize;
                loop {
                    let clen_offset = offset.min(self.current_len);
                    offset -= clen_offset;
                    self.current_len -= clen_offset;
                    self.current_ptr = unsafe { self.current_ptr.add(clen_offset) };
                    self.current_file_position += clen_offset;

                    if offset == 0 {
                        break Ok(self.current_file_position as u64);
                    }

                    if self.current_chunk_index >= self.chunks_count - 1 {
                        break Err(std::io::Error::new(
                            ErrorKind::UnexpectedEof,
                            "Unexpected eof",
                        ));
                    }

                    self.current_chunk_index += 1;
                    self.set_chunk_info(self.current_chunk_index);
                }
            }
            _ => {
                unimplemented!()
            }
        }
    }

    fn stream_position(&mut self) -> io::Result<u64> {
        let mut position = 0;

        let file_read = self.file.read();

        for i in 0..self.current_chunk_index {
            position += file_read.get_chunk(i).read().get_length();
        }

        position += file_read
            .get_chunk(self.current_chunk_index)
            .read()
            .get_length()
            - self.current_len;

        Ok(position as u64)
    }
}
