// Foreign-interface checks of the ThermoCore C ABI from the .NET Framework CLR
// (the runtime DWSIM uses on Windows). No DWSIM installation is needed.
// Usage: BridgeClientTests.exe <thermocore_bridge.dll> <request.json> <cli-molar-vapor-fraction>
using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Threading.Tasks;
using System.Web.Script.Serialization;
using ThermoCore.Bridge;

internal static class BridgeClientTests
{
    private static int failures;

    private static void Check(bool ok, string what)
    {
        Console.WriteLine((ok ? "PASS " : "FAIL ") + what);
        if (!ok) failures++;
    }

    private static string Expect(Action action, string kind)
    {
        try
        {
            action();
            return "no exception";
        }
        catch (ThermoCoreException e)
        {
            return e.Kind == kind ? null : "kind " + e.Kind + " (" + e.Message + ")";
        }
    }

    private static Dictionary<string, object> Load(string path)
    {
        return new JavaScriptSerializer().Deserialize<Dictionary<string, object>>(File.ReadAllText(path));
    }

    private static int Main(string[] args)
    {
        if (args.Length != 3)
        {
            Console.Error.WriteLine("usage: BridgeClientTests <thermocore_bridge.dll> <request.json> <cli-molar-vapor-fraction>");
            return 2;
        }
        ThermoCoreBridge.Load(args[0]);
        double cliBeta = double.Parse(args[2], CultureInfo.InvariantCulture);
        Check(ThermoCoreBridge.AbiVersion() == 1, "ABI version 1");

        var info = ThermoCoreBridge.Info();
        Check((string)info["schema"] == "thermocore.bridge.info.v1", "info schema");
        bool methaneCas = false;
        foreach (object c in (System.Collections.ArrayList)info["components"])
        {
            var d = ThermoCoreBridge.Obj(c);
            if ((string)d["id"] == "methane" && (string)d["cas"] == "74-82-8") methaneCas = true;
        }
        Check(methaneCas, "catalogue maps methane by CAS 74-82-8");

        var request = Load(args[1]);
        var r = ThermoCoreBridge.Calculate(request);
        double beta = ThermoCoreBridge.Num(r["molar_vapor_fraction"]);
        Check(Convert.ToInt32(r["phase_count"]) == 2, "two-phase split");
        Check(ThermoCoreBridge.Num(r["temperature_K"]) == 220.0 && ThermoCoreBridge.Num(r["pressure_Pa"]) == 1500000.0,
            "SI state echoed (K, Pa absolute)");
        // Same engine, same request: equal to the existing CLI interface up to
        // .NET Framework decimal parsing (about one ulp).
        Check(Math.Abs(beta - cliBeta) <= 4e-16, "vapor fraction equals thermo-cli: " + beta.ToString("R", CultureInfo.InvariantCulture) + " vs " + cliBeta.ToString("R", CultureInfo.InvariantCulture));
        var conv = ThermoCoreBridge.Obj(r["convergence"]);
        Check(ThermoCoreBridge.Num(conv["max_abs_ln_fugacity_residual"]) <= 1e-10, "fugacity residual <= 1e-10");
        var prov = ThermoCoreBridge.Obj(r["provenance"]);
        var model = ThermoCoreBridge.Obj(prov["model"]);
        Check(((string)model["formulation"]).Contains("zero-kij"), "interaction source survives in model formulation");
        Check(ThermoCoreBridge.Obj(prov["caloric_reference"])["definition"] is string, "caloric reference definition present");
        Check(ThermoCoreBridge.Obj(prov["data"])["version"] is string, "data version present");
        double sum = 0;
        foreach (object p in (System.Collections.ArrayList)r["phases"])
        {
            var ph = ThermoCoreBridge.Obj(p);
            sum += ThermoCoreBridge.Num(ph["molar_phase_fraction"]);
            Check(ThermoCoreBridge.Num(ph["mass_density_kg_per_m3"]) > 0, ph["label"] + " density in kg/m3 > 0");
        }
        Check(Math.Abs(sum - 1.0) < 1e-12, "phase fractions sum to one");

        // Typed failures surface as exceptions with the stable kind.
        string e;
        e = Expect(delegate
        {
            var q = Load(args[1]);
            ThermoCoreBridge.Obj(q["state"])["temperature_K"] = 150.0;
            ThermoCoreBridge.Calculate(q);
        }, "unsupported-capability");
        Check(e == null, "150 K outside package interval -> unsupported-capability " + e);
        e = Expect(delegate
        {
            var q = Load(args[1]);
            ((System.Collections.ArrayList)ThermoCoreBridge.Obj(q["package"])["components"])[4] = "hydrogen";
            ThermoCoreBridge.Calculate(q);
        }, "missing-parameter");
        Check(e == null, "unknown component -> missing-parameter " + e);
        e = Expect(delegate
        {
            var q = Load(args[1]);
            ThermoCoreBridge.Obj(q["state"])["pressure_bar"] = 15.0;
            ThermoCoreBridge.Calculate(q);
        }, "invalid-request");
        Check(e == null, "non-SI field rejected -> invalid-request " + e);
        string raw = ThermoCoreBridge.CalculateJson("{not json");
        Check(raw.Contains("\"status\":\"error\"") && raw.Contains("invalid-request"), "malformed JSON -> error document");
        bool rejected = false;
        try
        {
            var q = Load(args[1]);
            ThermoCoreBridge.Obj(q["state"])["pressure_Pa"] = double.NaN;
            ThermoCoreBridge.Calculate(q);
        }
        catch (ArgumentException) { rejected = true; }
        Check(rejected, "NaN rejected by the client before the call");

        // Reentrancy: concurrent calls give identical results.
        var results = new double[16];
        Parallel.For(0, results.Length, i =>
        {
            results[i] = ThermoCoreBridge.Num(ThermoCoreBridge.Calculate(Load(args[1]))["molar_vapor_fraction"]);
        });
        bool same = true;
        foreach (double b in results) same &= b == beta;
        Check(same, "16 concurrent calls identical");

        Console.WriteLine(failures == 0 ? "ALL PASSED" : failures + " FAILED");
        return failures == 0 ? 0 : 1;
    }
}
