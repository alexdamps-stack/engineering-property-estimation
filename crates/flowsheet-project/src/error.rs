use std::fmt;

/// Reasons a project cannot be read, written or turned into a package.
#[derive(Debug, Clone, PartialEq)]
pub enum ProjectError {
    /// Reading or writing the file failed.
    Io(String),
    /// The file is not valid JSON or does not match the project schema.
    Format(String),
    /// The file declares a `format_version` this program cannot read, or none.
    UnsupportedVersion(Option<u64>),
    /// The content is well formed but inconsistent (for example a feed
    /// composition for a unit that is not a feed).
    Invalid(String),
    /// The package selection cannot be built.
    Package(String),
}

impl fmt::Display for ProjectError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            ProjectError::Io(e) => write!(f, "file access failed: {e}"),
            ProjectError::Format(e) => write!(f, "not a valid project file: {e}"),
            ProjectError::UnsupportedVersion(Some(v)) => write!(
                f,
                "project format version {v} is not supported (this program reads version {})",
                crate::FORMAT_VERSION
            ),
            ProjectError::UnsupportedVersion(None) => {
                f.write_str("the file has no integer format_version")
            }
            ProjectError::Invalid(e) => write!(f, "inconsistent project: {e}"),
            ProjectError::Package(e) => write!(f, "property package: {e}"),
        }
    }
}

impl std::error::Error for ProjectError {}
