use std::fs::{self, File};
use std::io::{self, Write};
#[cfg(unix)]
use std::os::unix::fs::PermissionsExt;
use std::path::{Path, PathBuf};

use rsomics_common::{Context, Result, RsomicsError};
use tempfile::{Builder, NamedTempFile};

pub(crate) fn with_output_pair<T>(
    mapped: &Path,
    rejected: &Path,
    operation: impl FnOnce(&mut File, &mut File) -> Result<T>,
) -> Result<T> {
    let mut mapped_output = Staged::new(mapped)?;
    let mut rejected_output = Staged::new(rejected)?;
    let mapped_backup = Backup::new(mapped)?;
    let rejected_backup = Backup::new(rejected)?;
    let result = operation(
        mapped_output.temporary.as_file_mut(),
        rejected_output.temporary.as_file_mut(),
    )?;
    mapped_output.prepare()?;
    rejected_output.prepare()?;
    if let Err(error) = mapped_output.commit() {
        return Err(mapped_backup.restore(mapped, error));
    }
    if let Err(error) = rejected_output.commit() {
        let error = mapped_backup.restore(mapped, error);
        let error = rejected_backup.restore(rejected, error);
        return Err(error);
    }
    Ok(result)
}

enum Backup {
    Absent,
    Existing(NamedTempFile),
}

impl Backup {
    fn new(path: &Path) -> Result<Self> {
        match fs::metadata(path) {
            Ok(metadata) if !metadata.is_file() => Err(RsomicsError::InvalidInput(format!(
                "output {} is not a regular file",
                path.display()
            ))),
            Ok(_) => {
                let parent = path
                    .parent()
                    .filter(|parent| !parent.as_os_str().is_empty())
                    .unwrap_or_else(|| Path::new("."));
                let backup = Builder::new()
                    .prefix(".rsomics-liftover-backup-")
                    .make_in(parent, |backup| {
                        fs::hard_link(path, backup)?;
                        File::open(backup)
                    })
                    .rs_with_context(|| format!("backing up output {}", path.display()))?;
                Ok(Self::Existing(backup))
            }
            Err(error) if error.kind() == io::ErrorKind::NotFound => Ok(Self::Absent),
            Err(error) => Err(RsomicsError::Io(error)),
        }
    }

    fn restore(self, path: &Path, cause: RsomicsError) -> RsomicsError {
        let restored = match self {
            Self::Absent => match fs::remove_file(path) {
                Ok(()) => Ok(()),
                Err(error) if error.kind() == io::ErrorKind::NotFound => Ok(()),
                Err(error) => Err(error),
            },
            Self::Existing(backup) => backup
                .persist(path)
                .map(|_| ())
                .map_err(|error| error.error),
        };
        match restored {
            Ok(()) => cause,
            Err(error) => RsomicsError::Io(io::Error::new(
                error.kind(),
                format!(
                    "{cause}; also failed to restore output {}: {error}",
                    path.display()
                ),
            )),
        }
    }
}

struct Staged {
    path: PathBuf,
    parent: PathBuf,
    temporary: NamedTempFile,
}

impl Staged {
    fn new(path: &Path) -> Result<Self> {
        let parent = path
            .parent()
            .filter(|parent| !parent.as_os_str().is_empty())
            .unwrap_or_else(|| Path::new("."));
        let permissions = match fs::metadata(path) {
            Ok(metadata) => Some(metadata.permissions()),
            Err(error) if error.kind() == io::ErrorKind::NotFound => None,
            Err(error) => return Err(RsomicsError::Io(error)),
        };
        let mut builder = Builder::new();
        builder.prefix(".rsomics-liftover-");
        #[cfg(unix)]
        if permissions.is_none() {
            builder.permissions(fs::Permissions::from_mode(0o666));
        }
        if let Some(existing) = permissions.as_ref() {
            builder.permissions(existing.clone());
        }
        let temporary = builder
            .tempfile_in(parent)
            .rs_with_context(|| format!("creating output beside {}", path.display()))?;
        Ok(Self {
            path: path.to_owned(),
            parent: parent.to_owned(),
            temporary,
        })
    }

    fn prepare(&mut self) -> Result<()> {
        self.temporary
            .as_file_mut()
            .flush()
            .rs_with_context(|| format!("flushing output {}", self.path.display()))?;
        self.temporary
            .as_file_mut()
            .sync_all()
            .rs_with_context(|| format!("syncing output {}", self.path.display()))
    }

    fn commit(self) -> Result<()> {
        self.temporary.persist(&self.path).map_err(|error| {
            RsomicsError::Io(io::Error::new(
                error.error.kind(),
                format!("committing output {}: {}", self.path.display(), error.error),
            ))
        })?;
        #[cfg(unix)]
        File::open(&self.parent)
            .and_then(|directory| directory.sync_all())
            .rs_with_context(|| format!("syncing output directory {}", self.parent.display()))?;
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn second_commit_failure_rolls_back_first_output() {
        let directory = tempfile::tempdir().unwrap();
        let mapped = directory.path().join("mapped.bed");
        let rejected = directory.path().join("rejected.bed");
        fs::write(&mapped, b"old mapped\n").unwrap();
        fs::write(&rejected, b"old rejected\n").unwrap();

        let error = with_output_pair(&mapped, &rejected, |mapped_writer, rejected_writer| {
            mapped_writer.write_all(b"new mapped\n")?;
            rejected_writer.write_all(b"new rejected\n")?;
            fs::remove_file(&rejected)?;
            fs::create_dir(&rejected)?;
            Ok(())
        })
        .unwrap_err();

        assert!(error.to_string().contains("committing output"));
        assert_eq!(fs::read(&mapped).unwrap(), b"old mapped\n");
        assert!(rejected.is_dir());
    }
}
