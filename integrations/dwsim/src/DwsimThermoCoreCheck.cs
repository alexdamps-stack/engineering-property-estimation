// A12 DWSIM use case: independent ThermoCore equilibrium check of a DWSIM separator.
//
// A real DWSIM flowsheet (feed -> vessel -> vapor/liquid, DWSIM "Peng-Robinson (PR)")
// is built and solved headlessly through DWSIM's own Automation3 API. The feed state
// (K, Pa, mole fractions), the compound identities (by CAS number) and the PR
// interaction parameters DWSIM actually loaded are read back from the DWSIM objects
// and sent through the versioned ThermoCore bridge (C ABI v1). The ThermoCore result,
// with provenance and convergence evidence, is reported beside DWSIM's own split.
//
// DWSIM is a cross-check here, not a qualified reference: constants and the ideal-gas
// model differ. The consistency gates detect integration faults (units, ordering,
// kij mapping), not model accuracy.
//
// Usage: DwsimThermoCoreCheck.exe <dwsim-dir> <thermocore_bridge.dll> <report.json>
using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Reflection;
using System.Runtime.CompilerServices;
using System.Text;
using ThermoCore.Bridge;

internal static class DwsimThermoCoreCheck
{
    private static string dwsimDir;

    private static int Main(string[] args)
    {
        if (args.Length != 3)
        {
            Console.Error.WriteLine("usage: DwsimThermoCoreCheck <dwsim-dir> <thermocore_bridge.dll> <report.json>");
            return 2;
        }
        dwsimDir = Path.GetFullPath(args[0]);
        AppDomain.CurrentDomain.AssemblyResolve += delegate(object s, ResolveEventArgs e)
        {
            // DWSIM probes its optional ThermoC package from the ThermoCS subfolder.
            string name = new AssemblyName(e.Name).Name + ".dll";
            foreach (string dir in new[] { dwsimDir, Path.Combine(dwsimDir, "ThermoCS") })
            {
                string p = Path.Combine(dir, name);
                if (File.Exists(p)) return Assembly.LoadFrom(p);
            }
            return null;
        };
        Environment.SetEnvironmentVariable("PATH", dwsimDir + ";" + Environment.GetEnvironmentVariable("PATH"));
        string bridge = Path.GetFullPath(args[1]);
        string report = Path.GetFullPath(args[2]);
        Environment.CurrentDirectory = dwsimDir;
        try
        {
            return Run(bridge, report);
        }
        catch (Exception e)
        {
            Console.Error.WriteLine("FAILED: " + e);
            return 1;
        }
    }

    // Separate method so DWSIM types are resolved only after AssemblyResolve is attached.
    [MethodImpl(MethodImplOptions.NoInlining)]
    private static int Run(string bridgePath, string reportPath)
    {
        return Check.Execute(bridgePath, reportPath, dwsimDir);
    }
}

internal static class Check
{
    // DWSIM 9 compound names and the a priori consistency gates (absolute).
    private static readonly string[] Compounds = { "Methane", "Ethane", "Propane", "Nitrogen", "Carbon dioxide" };
    private static readonly double[] Feed = { 0.85, 0.075, 0.04666666666666667, 0.018333333333333333, 0.01 };
    private const double TemperatureK = 220.0, PressurePa = 1.5e6, MolarFlowMolPerS = 1.0;
    private const double GateVaporFraction = 0.02, GateComposition = 0.02, GateZ = 0.01;

    private static int failures;

    private static void Gate(bool ok, string what)
    {
        Console.WriteLine((ok ? "PASS " : "FAIL ") + what);
        if (!ok) failures++;
    }

    private static string F(double x)
    {
        return x.ToString("G9", CultureInfo.InvariantCulture);
    }

    public static int Execute(string bridgePath, string reportPath, string dwsimDir)
    {
        ThermoCoreBridge.Load(bridgePath);
        var info = ThermoCoreBridge.Info();
        var casToId = new Dictionary<string, string>();
        foreach (object c in (System.Collections.ArrayList)info["components"])
        {
            var d = ThermoCoreBridge.Obj(c);
            casToId[(string)d["cas"]] = (string)d["id"];
        }

        var automation = new DWSIM.Automation.Automation3();
        string dwsimVersion = automation.GetVersion();
        Console.WriteLine("DWSIM " + dwsimVersion + " from " + dwsimDir);
        DWSIM.Interfaces.IFlowsheet fs = automation.CreateFlowsheet();
        foreach (string n in Compounds) fs.AddCompound(n);

        string ppKey = null;
        foreach (string k in fs.AvailablePropertyPackages.Keys)
            if (k == "Peng-Robinson (PR)") ppKey = k;
        if (ppKey == null)
            throw new InvalidOperationException("DWSIM property package 'Peng-Robinson (PR)' not available; found: "
                + string.Join(", ", new List<string>(fs.AvailablePropertyPackages.Keys).ToArray()));
        DWSIM.Interfaces.IPropertyPackage pp = fs.CreateAndAddPropertyPackage(ppKey);

        var T = DWSIM.Interfaces.Enums.GraphicObjects.ObjectType.MaterialStream;
        var feed = (DWSIM.Interfaces.IMaterialStream)fs.AddObject(T, 50, 50, "FEED");
        var vap = (DWSIM.Interfaces.IMaterialStream)fs.AddObject(T, 250, 20, "VAPOR");
        var liq = (DWSIM.Interfaces.IMaterialStream)fs.AddObject(T, 250, 90, "LIQUID");
        var sep = fs.AddObject(DWSIM.Interfaces.Enums.GraphicObjects.ObjectType.Vessel, 150, 50, "SEPARATOR");
        fs.ConnectObjects(((DWSIM.Interfaces.ISimulationObject)feed).GraphicObject, sep.GraphicObject, 0, 0);
        fs.ConnectObjects(sep.GraphicObject, ((DWSIM.Interfaces.ISimulationObject)vap).GraphicObject, 0, 0);
        fs.ConnectObjects(sep.GraphicObject, ((DWSIM.Interfaces.ISimulationObject)liq).GraphicObject, 1, 0);
        foreach (var o in new DWSIM.Interfaces.ISimulationObject[] {
            (DWSIM.Interfaces.ISimulationObject)feed, (DWSIM.Interfaces.ISimulationObject)vap,
            (DWSIM.Interfaces.ISimulationObject)liq, sep })
            o.PropertyPackage = pp;

        feed.SetTemperature(TemperatureK);
        feed.SetPressure(PressurePa);
        feed.SetMolarFlow(MolarFlowMolPerS);
        feed.SetOverallComposition((double[])Feed.Clone());
        List<Exception> errors = automation.CalculateFlowsheet4(fs);
        if (errors != null && errors.Count > 0)
            throw new InvalidOperationException("DWSIM flowsheet calculation failed: " + errors[0]);

        // Read the state back from the solved DWSIM objects (DWSIM internal SI units).
        double tK = feed.GetTemperature(), pPa = feed.GetPressure();
        var feedPhase = feed.GetPhase("Mixture");
        var ids = new List<string>();
        var names = new List<string>();
        var z = new List<double>();
        foreach (var kv in feedPhase.Compounds)
        {
            string cas = kv.Value.ConstantProperties.CAS_Number;
            string id;
            if (!casToId.TryGetValue(cas, out id))
                throw new InvalidOperationException("DWSIM compound '" + kv.Key + "' (CAS " + cas + ") has no ThermoCore component");
            ids.Add(id);
            names.Add(kv.Key);
            z.Add(kv.Value.MoleFraction.Value);
        }
        double[][] kij = DwsimKij(pp, names);
        string source = dwsimVersion + " '" + ppKey + "' binary interaction parameters as loaded by the flowsheet property package (absent pairs are zero in DWSIM and passed as explicit zero)";

        var request = ThermoCoreBridge.TpFlashRequest("dwsim-separator-feed", "PR", ids, kij, source, z.ToArray(), tK, pPa);
        var tc = ThermoCoreBridge.Calculate(request);

        // DWSIM separator products.
        double nFeed = feed.GetMolarFlow(), nVap = vap.GetMolarFlow(), nLiq = liq.GetMolarFlow();
        double dwBeta = nVap / nFeed;
        double[] dwY = PhaseFractions(vap.GetPhase("Mixture"));
        double[] dwX = PhaseFractions(liq.GetPhase("Mixture"));
        double dwZv = vap.GetPhase("Vapor").Properties.compressibilityFactor.GetValueOrDefault(double.NaN);

        double tcBeta = ThermoCoreBridge.Num(tc["molar_vapor_fraction"]);
        var phases = (System.Collections.ArrayList)tc["phases"];
        var tcV = ThermoCoreBridge.Obj(phases[0]);
        var tcL = ThermoCoreBridge.Obj(phases[phases.Count - 1]);
        double[] tcY = ThermoCoreBridge.Nums(tcV["mole_fractions"]);
        double[] tcX = ThermoCoreBridge.Nums(tcL["mole_fractions"]);
        double tcZv = ThermoCoreBridge.Num(tcV["Z"]);
        double dy = MaxAbsDiff(dwY, tcY), dx = MaxAbsDiff(dwX, tcX);

        Console.WriteLine("feed  T=" + F(tK) + " K  P=" + F(pPa) + " Pa  components=" + string.Join(",", ids.ToArray()));
        Console.WriteLine("kij (from DWSIM):");
        foreach (double[] row in kij) Console.WriteLine("  [" + string.Join(", ", Array.ConvertAll(row, F)) + "]");
        Console.WriteLine("vapor fraction  DWSIM=" + F(dwBeta) + "  ThermoCore=" + F(tcBeta));
        Console.WriteLine("vapor Z         DWSIM=" + F(dwZv) + "  ThermoCore=" + F(tcZv));
        Console.WriteLine("max |dy|=" + F(dy) + "  max |dx|=" + F(dx));
        var model = ThermoCoreBridge.Obj(ThermoCoreBridge.Obj(tc["provenance"])["model"]);
        Console.WriteLine("ThermoCore model " + model["id"] + " " + model["implementation"]);

        Gate(Math.Abs(tK - TemperatureK) < 1e-9 && Math.Abs(pPa - PressurePa) < 1e-6, "DWSIM feed state read back in K and Pa");
        Gate(Math.Abs(nVap + nLiq - nFeed) <= 1e-9 * nFeed, "DWSIM separator mole balance closes");
        Gate(Convert.ToInt32(tc["phase_count"]) == 2 && dwBeta > 0 && dwBeta < 1, "both engines report a two-phase split");
        Gate(ThermoCoreBridge.Num(ThermoCoreBridge.Obj(tc["convergence"])["max_abs_ln_fugacity_residual"]) <= 1e-10, "ThermoCore fugacity residual <= 1e-10");
        Gate(((string)model["formulation"]).Contains(dwsimVersion), "DWSIM kij source survives into ThermoCore provenance");
        Gate(Math.Abs(dwBeta - tcBeta) <= GateVaporFraction, "|vapor fraction difference| <= " + F(GateVaporFraction));
        Gate(dy <= GateComposition && dx <= GateComposition, "phase composition differences <= " + F(GateComposition));
        Gate(Math.Abs(dwZv - tcZv) <= GateZ, "|vapor Z difference| <= " + F(GateZ));

        // Failure paths through the same host route.
        string kind = null;
        try
        {
            var bad = ThermoCoreBridge.TpFlashRequest("dwsim-cold", "PR", ids, kij, source, z.ToArray(), 150.0, pPa);
            ThermoCoreBridge.Calculate(bad);
        }
        catch (ThermoCoreException e) { kind = e.Kind; }
        Gate(kind == "unsupported-capability", "150 K state rejected with typed error (" + kind + ")");

        string unmapped = null;
        DWSIM.Interfaces.IFlowsheet fs2 = automation.CreateFlowsheet();
        fs2.AddCompound("Methane");
        fs2.AddCompound("N-butane");
        foreach (var c in fs2.SelectedCompounds.Values)
            if (!casToId.ContainsKey(c.CAS_Number)) unmapped = c.Name + " (CAS " + c.CAS_Number + ")";
        Gate(unmapped != null, "compound outside the ThermoCore catalogue detected before any call: " + unmapped);

        string flowsheetPath = Path.ChangeExtension(reportPath, ".dwxmz");
        automation.SaveFlowsheet(fs, flowsheetPath, true);

        var reportDoc = new Dictionary<string, object>
        {
            { "schema", "thermocore.a12.dwsim-separator-check.v1" },
            { "dwsim", new Dictionary<string, object> {
                { "version", dwsimVersion }, { "property_package", ppKey }, { "compounds", names },
                { "flowsheet_file", Path.GetFileName(flowsheetPath) },
                { "feed", new Dictionary<string, object> { { "temperature_K", tK }, { "pressure_Pa", pPa }, { "molar_flow_mol_per_s", nFeed }, { "mole_fractions", z } } },
                { "molar_vapor_fraction", dwBeta }, { "vapor_mole_fractions", dwY }, { "liquid_mole_fractions", dwX }, { "vapor_Z", dwZv } } },
            { "thermocore_request", request },
            { "thermocore_result", tc },
            { "deviation", new Dictionary<string, object> {
                { "molar_vapor_fraction", tcBeta - dwBeta }, { "max_abs_vapor_mole_fraction", dy },
                { "max_abs_liquid_mole_fraction", dx }, { "vapor_Z", tcZv - dwZv } } },
            { "gates", new Dictionary<string, object> {
                { "molar_vapor_fraction", GateVaporFraction }, { "mole_fraction", GateComposition }, { "Z", GateZ },
                { "purpose", "integration-fault detection (units, order, kij mapping); not reference qualification" } } },
            { "failed_gates", failures },
        };
        File.WriteAllText(reportPath, ThermoCoreBridge.ToJson(reportDoc), new UTF8Encoding(false));
        automation.ReleaseResources();
        Console.WriteLine(failures == 0 ? "ALL PASSED" : failures + " FAILED");
        return failures == 0 ? 0 : 1;
    }

    private static double[] PhaseFractions(DWSIM.Interfaces.IPhase phase)
    {
        var list = new List<double>();
        foreach (var kv in phase.Compounds) list.Add(kv.Value.MoleFraction.GetValueOrDefault(double.NaN));
        return list.ToArray();
    }

    private static double MaxAbsDiff(double[] a, double[] b)
    {
        if (a.Length != b.Length) return double.PositiveInfinity;
        double m = 0;
        for (int i = 0; i < a.Length; i++) m = Math.Max(m, Math.Abs(a[i] - b[i]));
        return m;
    }

    // Reads the kij matrix the DWSIM PR package uses (m_pr.InteractionParameters,
    // keyed by compound name in both orders). There is no public interface for it.
    private static double[][] DwsimKij(DWSIM.Interfaces.IPropertyPackage pp, IList<string> names)
    {
        const BindingFlags Any = BindingFlags.Public | BindingFlags.NonPublic | BindingFlags.Instance;
        FieldInfo f = pp.GetType().GetField("m_pr", Any);
        if (f == null) throw new InvalidOperationException("DWSIM PR package has no m_pr field; unsupported DWSIM version");
        object aux = f.GetValue(pp);
        var ip = (System.Collections.IDictionary)aux.GetType().GetProperty("InteractionParameters", Any).GetValue(aux, null);
        int n = names.Count;
        var k = new double[n][];
        for (int i = 0; i < n; i++)
        {
            k[i] = new double[n];
            for (int j = 0; j < n; j++)
                if (i != j) k[i][j] = Lookup(ip, names[i], names[j]);
        }
        for (int i = 0; i < n; i++)
            for (int j = 0; j < i; j++)
            {
                if (k[i][j] != k[j][i])
                {
                    if (k[i][j] != 0 && k[j][i] != 0)
                        throw new InvalidOperationException("DWSIM kij for " + names[i] + "/" + names[j] + " is asymmetric; ThermoCore requires a symmetric matrix");
                    // DWSIM stores one direction; use the one present.
                    double v = k[i][j] != 0 ? k[i][j] : k[j][i];
                    k[i][j] = v;
                    k[j][i] = v;
                }
            }
        return k;
    }

    private static double Lookup(System.Collections.IDictionary ip, string a, string b)
    {
        foreach (string key in new[] { a, a.ToLowerInvariant() })
        {
            if (!ip.Contains(key)) continue;
            var inner = (System.Collections.IDictionary)ip[key];
            foreach (string k2 in new[] { b, b.ToLowerInvariant() })
                if (inner.Contains(k2))
                {
                    object data = inner[k2];
                    return (double)data.GetType().GetProperty("kij").GetValue(data, null);
                }
        }
        return 0.0;
    }
}
