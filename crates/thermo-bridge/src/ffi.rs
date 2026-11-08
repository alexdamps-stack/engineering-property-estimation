//! C ABI of bridge v1; see `include/thermocore_bridge.h`.
//!
//! Ownership: every non-null `char*` returned by this library is a
//! NUL-terminated UTF-8 JSON document allocated here and must be released
//! exactly once with [`tc_bridge_free`]. Input pointers are borrowed for the
//! duration of the call. Functions are reentrant (no global mutable state).
//! A panic never unwinds into the host: it becomes an `internal-panic` error
//! response.
#![allow(unsafe_code)]

use crate::{ABI_VERSION, BridgeError, calculate, error_response, info};
use std::ffi::{CStr, CString, c_char};
use std::panic::{AssertUnwindSafe, catch_unwind};

fn into_c(value: serde_json::Value) -> *mut c_char {
    // serde_json escapes control characters, so the document has no interior NUL.
    CString::new(value.to_string())
        .map(CString::into_raw)
        .unwrap_or(std::ptr::null_mut())
}

fn guarded(f: impl FnOnce() -> serde_json::Value) -> *mut c_char {
    let value = catch_unwind(AssertUnwindSafe(f)).unwrap_or_else(|_| {
        error_response(
            None,
            &BridgeError {
                kind: "internal-panic",
                message: "the bridge caught an internal panic; no result was produced".into(),
            },
        )
    });
    into_c(value)
}

/// ABI version; hosts must refuse to proceed unless it equals the version they target.
#[unsafe(no_mangle)]
pub extern "C" fn tc_bridge_abi_version() -> u32 {
    ABI_VERSION
}

/// Bridge identity, capabilities and component catalogue (`thermocore.bridge.info.v1`).
#[unsafe(no_mangle)]
pub extern "C" fn tc_bridge_info() -> *mut c_char {
    guarded(info)
}

/// Evaluates one `thermocore.bridge.request.v1` document and returns a
/// `thermocore.bridge.response.v1` document. A null or non-UTF-8 request
/// yields an `invalid-request` response, never a null pointer.
///
/// # Safety
/// `request` must be null or point to a NUL-terminated byte string that stays
/// valid and unmodified for the duration of the call.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn tc_bridge_calculate(request: *const c_char) -> *mut c_char {
    guarded(|| {
        if request.is_null() {
            return error_response(None, &BridgeError::request("null request pointer"));
        }
        // SAFETY: non-null and NUL-terminated per the documented contract.
        let bytes = unsafe { CStr::from_ptr(request) };
        match bytes.to_str() {
            Ok(text) => calculate(text),
            Err(_) => error_response(None, &BridgeError::request("request is not valid UTF-8")),
        }
    })
}

/// Releases a document returned by this library. Null is ignored.
///
/// # Safety
/// `document` must be null or a pointer returned by `tc_bridge_info` or
/// `tc_bridge_calculate` that has not been freed yet.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn tc_bridge_free(document: *mut c_char) {
    if !document.is_null() {
        // SAFETY: allocated by CString::into_raw in this library, freed once.
        drop(unsafe { CString::from_raw(document) });
    }
}
