"""Reproduce GERG parameter extraction and independent NIST property fixtures.

Run with Python 3 and g++ on Linux. Never calls the Rust implementation.
Existing outputs require explicit --write; review diffs after regeneration.
"""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import tempfile

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "data/source-extracts/gerg"
REVISION = "3bdb9ab8ff317c618b0b59d1b704c2c86ddc5fce"
NAMES = ["methane", "nitrogen", "carbon-dioxide", "ethane", "propane",
         "isobutane", "n-butane", "isopentane", "n-pentane", "n-hexane",
         "n-heptane", "n-octane", "n-nonane", "n-decane", "hydrogen",
         "oxygen", "carbon-monoxide", "water", "hydrogen-sulfide", "helium", "argon"]
ARRAYS = "MMiGERG Dc Tc kpol kexp coik doik noik toik mNumb kpolij kexpij cijk eijk gijk nijk tijk dijk btij bvij gtij gvij fij th0i n0i".split()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    # Include the unmodified upstream translation so static parameter arrays
    # are visible after SetupGERG applies its documented transformations.
    preamble = '#include <iomanip>\n#include "GERG2008.cpp"\n'
    dumper = preamble + r'''
template<class T> void dump(const T& v) { std::cout << v; }
template<class T, size_t N> void dump(const T (&a)[N]) {
  std::cout << "[";
  for(size_t i=0;i<N;i++){if(i) std::cout << ","; dump(a[i]);}
  std::cout << "]";
}
int main(){SetupGERG(); std::cout << std::setprecision(17) << "{";
'''
    dumper += '\n'.join(f'std::cout << "{"," if i else ""}\\\"{name}\\\":"; dump({name});'
                        for i, name in enumerate(ARRAYS))
    dumper += '\nstd::cout << "}";}'
    runner = preamble + r'''
int main(){
  SetupGERG(); double T,D; std::cout << std::setprecision(17);
  while(std::cin >> T >> D){
    std::vector<double> x(22,0); for(int i=1;i<=21;i++) std::cin >> x[i];
    double p,z,pd,pdd,pdt,pt,u,h,s,cv,cp,w,g,jt,k,a;
    PropertiesGERG(T,D,x,p,z,pd,pdd,pdt,pt,u,h,s,cv,cp,w,g,jt,k,a);
    double ar[4][4],a0[3],tr,dr;
    AlpharGERG(1,0,T,D,x,ar); Alpha0GERG(T,D,x,a0);
    ReducingParametersGERG(x,tr,dr);
    double vals[]={p*1000,z,pd,pt*1000,u,h,s,cv,cp,w,a,g,
      ar[0][0],ar[0][1],ar[0][2],ar[1][0],ar[2][0],ar[1][1],a0[0],a0[1],a0[2],tr,dr*1000};
    for(double v:vals) std::cout << v << " "; std::cout << "\n";
  }
}
'''
    with tempfile.TemporaryDirectory(prefix="gerg-reference-") as folder:
        tmp = Path(folder)
        def compile_run(name, code, inputs=""):
            path = tmp / f"{name}.cpp"
            path.write_text(code)
            # SetupGERG transforms every array cell, including unused cells
            # of its local bijk table. Zero-init avoids indeterminate exports;
            # all coefficients actually used by the EOS are assigned upstream.
            subprocess.run(["g++", "-std=c++11", "-O2", "-ftrivial-auto-var-init=zero",
                            "-I", str(SOURCE), str(path), "-o", str(tmp/name)], check=True)
            return subprocess.check_output([str(tmp/name)], input=inputs, text=True)
        tables = json.loads(compile_run("tables", dumper))
        cases = []
        def add(name, x, t, rho):
            cases.append(dict(id=name, composition=dict(zip(NAMES, x)), temperature_k=t, molar_density=rho))
        # Every pure fluid, all departure-function families, and all 210 pairs.
        for i, name in enumerate(NAMES):
            x = [float(j == i) for j in range(21)]
            add(f"pure-{name}", x, 400., 1.)
            add(f"dense-{name}", x, 450., 1000.)
        for i in range(21):
            for j in range(i+1,21):
                x = [0.]*21
                x[i], x[j] = .7, .3
                add(f"pair-{NAMES[i]}-{NAMES[j]}", x, 400., 100.)
        natural = [.9,.02,.02,.04,.02]+[0.]*16
        for t in [200.,250.,300.,350.,450.]:
            for rho in [1.,100.,1000.,4000.]:
                add(f"natural-gas-{t}-{rho}", natural, t, rho)
        all21 = [.77824,.02,.06,.08,.03,.0015,.003,.0005,.00165,.00215,.00088,.00024,.00015,.00009,.004,.005,.002,.0001,.0025,.007,.001]
        add("nist-all21", all21, 400., 1000.)
        inputs = "\n".join(" ".join(map(str,[c["temperature_k"],c["molar_density"]/1000,*c["composition"].values()])) for c in cases)
        keys = "pressure_pa compressibility dp_drho dp_dt internal_energy enthalpy entropy cv cp speed_of_sound helmholtz_energy gibbs_energy ar ar_d ar_dd ar_t ar_tt ar_dt a0 a0_t a0_tt reducing_temperature reducing_density".split()
        for c, row in zip(cases, compile_run("reference", runner, inputs).splitlines(), strict=True):
            c["expected"] = dict(zip(keys,map(float,row.split()),strict=True))
            c["composition"] = {k:v for k,v in c["composition"].items() if v}
        provenance = dict(upstream="https://github.com/usnistgov/AGA8", revision=REVISION,
                          source_sha256=hashlib.sha256((SOURCE/"GERG2008.cpp").read_bytes()).hexdigest(),
                          generator_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
        outputs = {
            ROOT/"data/gerg-2008.json": dict(provenance=provenance, components=NAMES, tables=tables),
            ROOT/"validation/fixtures/gerg-v1.json": dict(provenance=provenance,
                reference="NIST AGA8 GERG2008 C++ 2.01; native ideal reference; R=8.314472 J/(mol K)",
                tolerance=dict(relative=2e-10, absolute=2e-8), cases=cases),
        }
        for path, content in outputs.items():
            encoded = json.dumps(content, indent=2, allow_nan=False)+"\n"
            if args.write:
                path.write_text(encoded)
            else:
                assert path.read_text() == encoded, f"reproduction differs: {path}"
        print(f"{'Wrote' if args.write else 'Verified'} parameters and {len(cases)} independent states")


if __name__ == "__main__":
    main()
