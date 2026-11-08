"""Independent GERG pressure/density references: NIST density + teqp/SciPy.

By default verifies the committed fixture without writing. An explicit --output
writes a candidate. Never invokes Rust. Uses reference-requirements.txt and g++.
"""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import tempfile
import numpy as np
import scipy
from scipy.optimize import brentq
import teqp

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT/"data/source-extracts/gerg"
NAMES = json.loads((ROOT/"data/gerg-2008.json").read_text())["components"]
ALIASES = {"carbon-dioxide":"carbondioxide", "carbon-monoxide":"carbonmonoxide", "hydrogen-sulfide":"hydrogensulfide"}
assert teqp.__version__ == "0.23.2"
R = 8.314472


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    cases = []
    def add(name, composition, t, p, guess=0.):
        cases.append(dict(id=name, composition=composition, temperature_k=t, pressure_pa=p, nist_initial_density=guess))
    for name in NAMES:
        add("pure-"+name, {name:1.},400.,100000.)
    gas = dict(methane=.9, ethane=.04, propane=.02, nitrogen=.02, **{"carbon-dioxide":.02})
    for t in [250.,300.,350.,450.]:
        for p in [100000.,1e6,1e7,3e7]:
            add(f"natural-gas-{t}-{p}", gas,t,p)
    # Dense propane at this low pressure needs a close positive-pressure seed
    # for NIST's logarithmic-pressure Newton routine. Independently located by
    # teqp; root agreement is checked again with a fresh SciPy bracket below.
    for name,t,p,liquid in [("propane",250.,1e5,12664.),("ethane",250.,1e6,14850.)]:
        add(name+"-low-density",{name:1.},t,p)
        add(name+"-high-density",{name:1.},t,p,-liquid)
    # Full GERG component set, including new GERG-2008 species.
    full = [.77824,.02,.06,.08,.03,.0015,.003,.0005,.00165,.00215,.00088,.00024,.00015,.00009,.004,.005,.002,.0001,.0025,.007,.001]
    add("all21-400K-20MPa",dict(zip(NAMES,full)),400.,2e7)
    code = r'''
#include <iomanip>
#include "GERG2008.cpp"
int main(){SetupGERG(); double T,P,D; int row=0;
  std::cout << std::setprecision(17);
  while(std::cin >> T >> P >> D){
    ++row;
    std::vector<double> x(22,0); for(int i=1;i<=21;i++) std::cin >> x[i];
    int ierr; std::string error;
    double initialD=D;
    DensityGERG(0,T,P/1000,x,D,ierr,error);
    if(ierr){double pp,zz;PressureGERG(T,std::abs(initialD),x,pp,zz);std::cerr << "case " << row << " T=" << T << " P=" << P << " Dguess=" << initialD << " pguess=" << pp << " stiffness=" << dPdDsave << ": " << error; return 1;}
    double pp,z,pd,pdd,pdt,pt,u,h,s,cv,cp,w,g,jt,k,a;
    PropertiesGERG(T,D,x,pp,z,pd,pdd,pdt,pt,u,h,s,cv,cp,w,g,jt,k,a);
    double vals[]={D*1000,pp*1000,z,h,s,cp,cv,w,pd,pt*1000};
    for(double v:vals) std::cout << v << " "; std::cout << "\n";
  }
}
'''
    inputs = "\n".join(" ".join(map(str,[c["temperature_k"],c["pressure_pa"],c["nist_initial_density"]/1000,
                            *[c["composition"].get(n,0.) for n in NAMES]])) for c in cases)
    with tempfile.TemporaryDirectory(prefix="gerg-density-reference-") as folder:
        tmp = Path(folder)
        (tmp/"reference.cpp").write_text(code)
        subprocess.run(["g++","-std=c++11","-O2","-ftrivial-auto-var-init=zero","-I",str(SOURCE),str(tmp/"reference.cpp"),"-o",str(tmp/"reference")],check=True)
        lines = subprocess.check_output([str(tmp/"reference")],input=inputs,text=True).splitlines()
    for c,line in zip(cases,lines,strict=True):
        expected = dict(zip("molar_density pressure_pa compressibility enthalpy entropy cp cv speed_of_sound dp_drho dp_dt".split(),map(float,line.split()),strict=True))
        rho = expected["molar_density"]
        assert 1e-3 < rho < 30000 and expected["dp_drho"] > 0 and expected["cv"] > 0, c
        assert abs(expected["pressure_pa"]-c["pressure_pa"]) <= 1e-7+2e-10*c["pressure_pa"], c
        # NIST supplies an independent density root; teqp AD pressure and SciPy
        # then independently solve a narrow bracket around that particular root.
        lo,hi = .96*rho,1.02*rho
        x = np.array(list(c["composition"].values()))
        model = teqp.make_model(dict(kind="GERG2008resid",model=dict(names=[ALIASES.get(n,n) for n in c["composition"]])))
        def pressure_residual(d):
            return d*R*c["temperature_k"]*(1+model.get_Ar01(c["temperature_k"],d,x))-c["pressure_pa"]
        independent = brentq(pressure_residual,lo,hi,xtol=1e-10,rtol=2e-14)
        assert abs(independent-rho) <= 2e-7+2e-9*rho, (c["id"],independent,rho)
        c.update(density_bracket=[lo,hi],expected=expected,teqp_molar_density=independent)
    output = dict(reference="NIST AGA8 C++ GERG2008 2.01 density solver, independent teqp 0.23.2 AD/SciPy Brent roots",
        versions=dict(teqp=teqp.__version__,numpy=np.__version__,scipy=scipy.__version__),
        source_revision="3bdb9ab8ff317c618b0b59d1b704c2c86ddc5fce",
        source_sha256=hashlib.sha256((SOURCE/"GERG2008.cpp").read_bytes()).hexdigest(),
        generator_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        tolerances={"molar_density":[2e-7,2e-9],"pressure_pa":[1e-5,2e-9],"compressibility":[2e-9,2e-9],
            "enthalpy":[1e-5,2e-9],"entropy":[2e-7,2e-9],"cp":[2e-6,2e-8],"cv":[2e-6,2e-8],
            "speed_of_sound":[2e-4,2e-8],"dp_drho":[2e-5,2e-8],"dp_dt":[2e-5,2e-8]},cases=cases)
    if args.output:
        args.output.write_text(json.dumps(output,indent=2,allow_nan=False)+"\n")
        print(f"Wrote {len(cases)} independently solved GERG TP roots to {args.output}")
    else:
        assert json.loads((ROOT/"validation/fixtures/gerg-density-v1.json").read_text()) == output, "GERG density reproduction differs; review a candidate explicitly"
        print(f"Verified {len(cases)} independently solved GERG TP roots")


if __name__ == "__main__":
    main()
