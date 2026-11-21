//! C ABI of bridge v1; see `include/thermocore_bridge.h`.
//!
//! Ownership: every non-null `char*` returned by this library is a
//! NUL-terminated UTF-8 JSON document allocated here and must be released
//! exactly once with [`tc_bridge_free`]. Input pointers are borrowed for the
//! duration of the call. Functions are reentrant (no global mutable state).
//! A panic never unwinds into the host: it becomes an `internal-panic` error
//! response, which echoes the request's `request_id` when the request was
//! parsed before the panic.
#![allow(unsafe_code)]

use crate::{ABI_VERSION, BridgeError, calculate_recording_id, error_response, info};
use std::ffi::{CStr, CString, c_char};
use std::panic::{AssertUnwindSafe, catch_unwind};

fn into_c(value: serde_json::Value) -> *mut c_char {
    // serde_json escapes control characters, so the document has no interior NUL.
    CString::new(value.to_string())
        .map(CString::into_raw)
        .unwrap_or(std::ptr::null_mut())
}

/// Runs `f`, turning a panic into an `internal-panic` response. `f` may store
/// the request's `request_id` in the slot it is given; the panic response
/// echoes whatever was stored there before the panic.
fn guarded(f: impl FnOnce(&mut Option<String>) -> serde_json::Value) -> *mut c_char {
    into_c(guarded_value(f))
}

fn guarded_value(f: impl FnOnce(&mut Option<String>) -> serde_json::Value) -> serde_json::Value {
    let mut request_id = None;
    let outcome = catch_unwind(AssertUnwindSafe(|| f(&mut request_id)));
    outcome.unwrap_or_else(|_| {
        error_response(
            request_id.as_deref(),
            &BridgeError {
                kind: "internal-panic",
                message: "the bridge caught an internal panic; no result was produced".into(),
            },
        )
    })
}

/// ABI version; hosts must refuse to proceed unless it equals the version they target.
#[unsafe(no_mangle)]
pub extern "C" fn tc_bridge_abi_version() -> u32 {
    ABI_VERSION
}

/// Bridge identity, capabilities and component catalogue (`thermocore.bridge.info.v1`).
#[unsafe(no_mangle)]
pub extern "C" fn tc_bridge_info() -> *mut c_char {
    guarded(|_| info())
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
    guarded(|request_id| {
        if request.is_null() {
            return error_response(None, &BridgeError::request("null request pointer"));
        }
        // SAFETY: non-null and NUL-terminated per the documented contract.
        let bytes = unsafe { CStr::from_ptr(request) };
        match bytes.to_str() {
            Ok(text) => calculate_recording_id(text, request_id),
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

#[cfg(test)]
mod tests {
    use super::guarded_value;

    #[test]
    fn panic_response_keeps_a_recorded_request_id() {
        let r = guarded_value(|id| {
            *id = Some("req-7".into());
            panic!("simulated evaluation panic");
        });
        assert_eq!(r["status"], "error");
        assert_eq!(r["error"]["kind"], "internal-panic");
        assert_eq!(r["request_id"], "req-7");
    }

    #[test]
    fn panic_before_parsing_has_null_request_id() {
        let r = guarded_value(|_| panic!("simulated early panic"));
        assert_eq!(r["error"]["kind"], "internal-panic");
        assert!(r["request_id"].is_null());
    }

    #[test]
    fn calculate_records_the_request_id_before_evaluating() {
        // Parses as a v1 document but fails in evaluation (wrong schema), so
        // the id must already be recorded by then.
        let mut id = None;
        let r = crate::calculate_recording_id(
            r#"{"schema":"other","request_id":"r1","calculation":"tp-flash",
                "package":{"model":"PR","components":[],"kij":[],"interaction_source":"x"},
                "state":{"mole_fractions":[]}}"#,
            &mut id,
        );
        assert_eq!(r["status"], "error");
        assert_eq!(id.as_deref(), Some("r1"));
    }
}
