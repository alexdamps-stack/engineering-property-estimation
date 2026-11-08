/* ThermoCore bridge, C ABI version 1 (schema thermocore.bridge.*.v1).
 * See docs/a12-status.md for the JSON contract. All strings are UTF-8 JSON.
 * Every returned pointer is owned by the library: release it with
 * tc_bridge_free exactly once. Functions are reentrant. */
#ifndef THERMOCORE_BRIDGE_H
#define THERMOCORE_BRIDGE_H
#include <stdint.h>
#ifdef __cplusplus
extern "C" {
#endif

#define TC_BRIDGE_ABI_VERSION 1u

uint32_t tc_bridge_abi_version(void);
char *tc_bridge_info(void);
char *tc_bridge_calculate(const char *request_json);
void tc_bridge_free(char *document);

#ifdef __cplusplus
}
#endif
#endif
