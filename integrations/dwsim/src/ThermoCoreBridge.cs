// ThermoCore bridge client for .NET Framework 4.6.2+ hosts (DWSIM on Windows).
// Wraps the C ABI v1 of thermocore_bridge.dll (crates/thermo-bridge). The
// library has no DWSIM dependency; DwsimThermoCoreCheck.cs adapts DWSIM objects.
// Written for the C# 5 compiler shipped with .NET Framework (csc.exe v4.0.30319).
using System;
using System.Collections;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Runtime.InteropServices;
using System.Text;
using System.Web.Script.Serialization;

namespace ThermoCore.Bridge
{
    /// <summary>A typed bridge error; Kind is the stable v1 error category.</summary>
    public sealed class ThermoCoreException : Exception
    {
        public string Kind { get; private set; }
        public string RequestId { get; private set; }

        public ThermoCoreException(string kind, string message, string requestId)
            : base("ThermoCore " + kind + ": " + message)
        {
            Kind = kind;
            RequestId = requestId;
        }
    }

    public static class ThermoCoreBridge
    {
        public const uint TargetAbiVersion = 1;
        public const string RequestSchema = "thermocore.bridge.request.v1";
        public const string ResponseSchema = "thermocore.bridge.response.v1";
        private const string Library = "thermocore_bridge";

        [DllImport(Library, CallingConvention = CallingConvention.Cdecl)]
        private static extern uint tc_bridge_abi_version();

        [DllImport(Library, CallingConvention = CallingConvention.Cdecl)]
        private static extern IntPtr tc_bridge_info();

        [DllImport(Library, CallingConvention = CallingConvention.Cdecl)]
        private static extern IntPtr tc_bridge_calculate(byte[] requestUtf8);

        [DllImport(Library, CallingConvention = CallingConvention.Cdecl)]
        private static extern void tc_bridge_free(IntPtr document);

        [DllImport("kernel32", SetLastError = true, CharSet = CharSet.Unicode)]
        private static extern IntPtr LoadLibrary(string path);

        private static readonly JavaScriptSerializer Json = new JavaScriptSerializer { MaxJsonLength = int.MaxValue };
        private static bool loaded;

        /// <summary>
        /// Loads thermocore_bridge.dll from an explicit path (so the host's DLL
        /// search order is irrelevant) and refuses any other ABI version.
        /// </summary>
        public static void Load(string dllPath)
        {
            string full = Path.GetFullPath(dllPath);
            if (!File.Exists(full))
                throw new FileNotFoundException("ThermoCore bridge library not found", full);
            if (LoadLibrary(full) == IntPtr.Zero)
                throw new InvalidOperationException("LoadLibrary failed for " + full + " (Win32 error " + Marshal.GetLastWin32Error() + ")");
            uint abi = tc_bridge_abi_version();
            if (abi != TargetAbiVersion)
                throw new InvalidOperationException("ThermoCore bridge ABI " + abi + " is incompatible with client ABI " + TargetAbiVersion);
            loaded = true;
        }

        public static uint AbiVersion()
        {
            return tc_bridge_abi_version();
        }

        private static string Take(IntPtr document)
        {
            if (document == IntPtr.Zero)
                throw new InvalidOperationException("ThermoCore bridge returned a null document");
            try
            {
                int n = 0;
                while (Marshal.ReadByte(document, n) != 0) n++;
                byte[] bytes = new byte[n];
                Marshal.Copy(document, bytes, 0, n);
                return Encoding.UTF8.GetString(bytes);
            }
            finally
            {
                tc_bridge_free(document);
            }
        }

        private static void EnsureLoaded()
        {
            if (!loaded)
                throw new InvalidOperationException("call ThermoCoreBridge.Load(path) first");
        }

        /// <summary>Raw info document (thermocore.bridge.info.v1).</summary>
        public static string InfoJson()
        {
            EnsureLoaded();
            return Take(tc_bridge_info());
        }

        public static Dictionary<string, object> Info()
        {
            return Json.Deserialize<Dictionary<string, object>>(InfoJson());
        }

        /// <summary>Raw v1 request/response. Never throws for calculation failures.</summary>
        public static string CalculateJson(string requestJson)
        {
            EnsureLoaded();
            byte[] utf8 = Encoding.UTF8.GetBytes(requestJson + "\0");
            return Take(tc_bridge_calculate(utf8));
        }

        /// <summary>Returns the v1 "result" object or throws ThermoCoreException.</summary>
        public static Dictionary<string, object> Calculate(IDictionary<string, object> request)
        {
            string text = CalculateJson(ToJson(request));
            var response = Json.Deserialize<Dictionary<string, object>>(text);
            if ((string)response["schema"] != ResponseSchema)
                throw new InvalidOperationException("unexpected response schema " + response["schema"]);
            string id = response["request_id"] as string;
            if ((string)response["status"] == "ok")
                return (Dictionary<string, object>)response["result"];
            var error = (Dictionary<string, object>)response["error"];
            throw new ThermoCoreException((string)error["kind"], (string)error["message"], id);
        }

        /// <summary>Builds a v1 TP-flash request; all quantities SI (K, Pa absolute).</summary>
        public static Dictionary<string, object> TpFlashRequest(string requestId, string model,
            IList<string> components, double[][] kij, string interactionSource,
            double[] moleFractions, double temperatureK, double pressurePa)
        {
            return new Dictionary<string, object>
            {
                { "schema", RequestSchema },
                { "request_id", requestId },
                { "calculation", "tp-flash" },
                { "package", new Dictionary<string, object>
                    {
                        { "model", model },
                        { "components", components },
                        { "kij", kij },
                        { "interaction_source", interactionSource },
                    } },
                { "state", new Dictionary<string, object>
                    {
                        { "mole_fractions", moleFractions },
                        { "temperature_K", temperatureK },
                        { "pressure_Pa", pressurePa },
                    } },
            };
        }

        /// <summary>
        /// Minimal JSON writer for requests. Doubles use "G17", which always
        /// round-trips; JavaScriptSerializer uses "R", which on .NET Framework
        /// x64 does not round-trip for some values. Non-finite numbers are
        /// rejected here rather than sent as invalid JSON.
        /// </summary>
        public static string ToJson(object value)
        {
            var sb = new StringBuilder();
            Write(sb, value);
            return sb.ToString();
        }

        private static void Write(StringBuilder sb, object v)
        {
            if (v == null) { sb.Append("null"); return; }
            string s = v as string;
            if (s != null) { WriteString(sb, s); return; }
            if (v is bool) { sb.Append((bool)v ? "true" : "false"); return; }
            if (v is double || v is float)
            {
                double d = Convert.ToDouble(v, CultureInfo.InvariantCulture);
                if (double.IsNaN(d) || double.IsInfinity(d))
                    throw new ArgumentException("non-finite number cannot be sent to ThermoCore");
                sb.Append(d.ToString("G17", CultureInfo.InvariantCulture));
                return;
            }
            // JavaScriptSerializer reads JSON numbers as decimal: keep their exact text.
            if (v is decimal) { sb.Append(((decimal)v).ToString(CultureInfo.InvariantCulture)); return; }
            if (v is int || v is long || v is uint) { sb.Append(Convert.ToString(v, CultureInfo.InvariantCulture)); return; }
            var dict = v as IDictionary;
            if (dict != null)
            {
                sb.Append('{');
                bool first = true;
                foreach (DictionaryEntry e in dict)
                {
                    if (!first) sb.Append(',');
                    first = false;
                    WriteString(sb, (string)e.Key);
                    sb.Append(':');
                    Write(sb, e.Value);
                }
                sb.Append('}');
                return;
            }
            var seq = v as IEnumerable;
            if (seq != null)
            {
                sb.Append('[');
                bool first = true;
                foreach (object x in seq)
                {
                    if (!first) sb.Append(',');
                    first = false;
                    Write(sb, x);
                }
                sb.Append(']');
                return;
            }
            throw new ArgumentException("unsupported JSON value type " + v.GetType());
        }

        private static void WriteString(StringBuilder sb, string s)
        {
            sb.Append('"');
            foreach (char c in s)
            {
                if (c == '"') sb.Append("\\\"");
                else if (c == '\\') sb.Append("\\\\");
                else if (c < 0x20) sb.Append("\\u").Append(((int)c).ToString("x4"));
                else sb.Append(c);
            }
            sb.Append('"');
        }

        // Small typed accessors for the untyped JavaScriptSerializer tree.
        // Note: .NET Framework number parsing is not guaranteed to be correctly
        // rounded, so values read here can differ from the Rust result by about
        // one unit in the last place; the Rust side emits shortest round-trip text.
        public static double Num(object o)
        {
            return Convert.ToDouble(o, CultureInfo.InvariantCulture);
        }

        public static double[] Nums(object o)
        {
            var list = new List<double>();
            foreach (object x in (IEnumerable)o) list.Add(Num(x));
            return list.ToArray();
        }

        public static Dictionary<string, object> Obj(object o)
        {
            return (Dictionary<string, object>)o;
        }
    }
}
