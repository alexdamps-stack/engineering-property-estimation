"""Independent GERG TP/PH references, never invoking Rust.

NIST C++ pressure/ideal calorics, teqp AD chemical potentials and SciPy
multi-start TPD minimization / Powell-hybrid phase splitting. Writes only an
explicit --output candidate. Reuses the independently implemented A4 algorithms.
"""
import argparse
import ctypes
import hashlib
import importlib.util
import importlib.metadata
import json
from pathlib import Path
import subprocess
import tempfile
import numpy as np
import teqp
from scipy.optimize import minimize, least_squares, root, brentq

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT/"data/source-extracts/gerg"
spec = importlib.util.spec_from_file_location("a4_reference",ROOT/"validation/generate-flash-references.py")
a4 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(a4)
a4.R = R = 8.314472
a4.RANDOM_STARTS = 16
a4.SEED = 0xA6A6
NAMES = json.loads((ROOT/"data/gerg-2008.json").read_text())["components"]
ALIASES = {"carbon-dioxide":"carbondioxide", "carbon-monoxide":"carbonmonoxide", "hydrogen-sulfide":"hydrogensulfide"}

CODE = r'''
#include <algorithm>
#include "GERG2008.cpp"
extern "C" void setup(){SetupGERG();}
void reset(){Told=0;Trold2=0;for(int i=1;i<=21;i++)xold[i]=-1;}
extern "C" int density_roots(double T,double p,const double* fractions,double* out){
  std::vector<double> x(22,0);for(int i=1;i<=21;i++)x[i]=fractions[i-1];reset();
  std::vector<double> grid;
  for(int i=0;i<=320;i++){double q=i/320.;grid.push_back(.001*exp(q*log(30000/.001)));grid.push_back(.001+q*(30000-.001));}
  std::sort(grid.begin(),grid.end());int count=0;double previous=grid[0],pp,z;
  double first_negative=1e100,last_negative=-1;
  PressureGERG(T,previous/1000,x,pp,z);double left=pp*1000-p;
  for(size_t j=1;j<grid.size();j++){
    double current=grid[j];PressureGERG(T,current/1000,x,pp,z);double right=pp*1000-p;
    if(dPdDsave<=0){first_negative=std::min(first_negative,current);last_negative=current;}
    if(left<=0 && right>=0){
      double lo=previous,hi=current,rho=lo;
      for(int it=0;it<100;it++){
        rho=(lo+hi)/2;PressureGERG(T,rho/1000,x,pp,z);double f=pp*1000-p;
        if(std::abs(f)<=2e-13*std::max(p,1.))break;
        if(f<0)lo=rho;else hi=rho;
      }
      double pd,pdd,pdt,pt,u,h,s,cv,cp,w,g,jt,k,a;
      PropertiesGERG(T,rho/1000,x,pp,z,pd,pdd,pdt,pt,u,h,s,cv,cp,w,g,jt,k,a);
      if(pd>0 && cv>0 && (count==0 || std::abs(rho-out[count-1])>1e-7*rho)){
        if(count>=16)return -1;out[count++]=rho;
      }
    }
    previous=current;left=right;
  }
  int outer=0;
  for(int i=0;i<count;i++)if(last_negative<0 || out[i]<=first_negative || out[i]>=last_negative)out[outer++]=out[i];
  return outer;
}
extern "C" void properties(double T,double rho,const double* fractions,double* out){
  std::vector<double> x(22,0);for(int i=1;i<=21;i++)x[i]=fractions[i-1];reset();
  double p,z,pd,pdd,pdt,pt,u,h,s,cv,cp,w,g,jt,k,a;
  PropertiesGERG(T,rho/1000,x,p,z,pd,pdd,pdt,pt,u,h,s,cv,cp,w,g,jt,k,a);
  out[0]=h;out[1]=s;out[2]=cp;out[3]=rho;out[4]=z;
}
'''


class Reference:
    def __init__(self,lib,ids):
        self.lib,self.ids=lib,ids
        self.model=teqp.make_model(dict(kind="GERG2008resid",model=dict(names=[ALIASES.get(n,n) for n in ids])))
    def canonical(self,x):
        return np.array([x[self.ids.index(n)] if n in self.ids else 0. for n in NAMES],dtype=float)
    def get_fugacity_coefficients(self,t,rhovec):
        return self.model.get_fugacity_coefficients(t,rhovec)
    def roots(self,t,p,x):
        canonical=self.canonical(x);out=np.zeros(16)
        count=self.lib.density_roots(t,p,canonical.ctypes.data_as(PTR),out.ctypes.data_as(PTR))
        assert count>0,(self.ids,t,p,x,"no regular density root")
        # Only outer roots: the interior two-phase region of multiparameter
        # pure EOS can contain nonphysical extra roots. No completeness proof.
        densities=[out[0]] if count==1 else [out[0],out[count-1]]
        result=[]
        for rho in densities:
            check=rho*R*t*(1+self.model.get_Ar01(t,rho,x))
            assert abs(check-p)<=2e-10*p,(check,p)
            phi=np.log(self.model.get_fugacity_coefficients(t,rho*x))
            result.append((p/(rho*R*t),rho,phi))
        return result
    def properties(self,t,rho,x):
        canonical=self.canonical(x);out=np.zeros(5)
        self.lib.properties(t,rho,canonical.ctypes.data_as(PTR),out.ctypes.data_as(PTR))
        return dict(zip(["h","s","cp","rho","Z"],out.tolist()))


PTR=ctypes.POINTER(ctypes.c_double)


def tpd_search(model,ids,t,p,z):
    """Independent constrained numerical evidence over regular exterior roots.

    Failed minimization trajectories are retained, never treated as stable
    points. Successful trajectories supply finite evidence, not global proof.
    """
    active=z>0
    feed=a4.min_gibbs(model,t,p,z)
    if active.sum()<2:return 0.,z.tolist(),feed[0]
    d=np.where(active,np.log(np.where(active,z,1.))+feed[2],0.)
    k=a4.wilson(ids,t,p)
    rng=np.random.default_rng(a4.SEED)
    starts=[z*k,z/k]+[rng.dirichlet(np.ones(len(z)))*active for _ in range(a4.RANDOM_STARTS)]
    def objective(logw):
        wbig=np.zeros(len(z));wbig[active]=np.exp(logw);w=wbig/wbig.sum()
        phi=a4.min_gibbs(model,t,p,w)[2]
        g=logw+phi[active]-d[active]
        return 1.+float(np.sum(wbig[active]*(g-1.))),wbig[active]*g
    best=(float("inf"),None)
    failures=[]
    for i,start in enumerate(starts):
        try:
            result=minimize(objective,np.log(np.maximum(start[active],1e-12)),jac=True,
                method="L-BFGS-B",bounds=[(-60.,5.)]*active.sum(),
                options=dict(ftol=1e-15,gtol=1e-12,maxiter=2000))
            assert np.isfinite(result.fun)
            if result.fun<best[0]:best=(float(result.fun),result.x)
            if not result.success:failures.append(dict(start=i,reason=str(result.message)))
        except AssertionError as error:
            failures.append(dict(start=i,reason=str(error)))
    assert best[1] is not None,"all reference TPD trajectories failed"
    model.trial_records.append(dict(T=t,p=p,composition=z.tolist(),attempted_starts=len(starts),failed_trajectories=failures))
    w=np.zeros(len(z));w[active]=np.exp(best[1]);w/=w.sum()
    return best[0],w.tolist(),a4.min_gibbs(model,t,p,w)[0]


def split(model,t,p,z,lnk0):
    """Independent trust-region least squares on ln K, with brentq RR."""
    active=z>0
    def residual(lnk):
        k=np.ones(len(z));k[active]=np.exp(lnk)
        _,x,y=a4.rachford_rice(z,k)
        return lnk-(a4.min_gibbs(model,t,p,x)[2]-a4.min_gibbs(model,t,p,y)[2])[active]
    failures=[]
    starts=[lnk0[active],np.log(a4.wilson(model.ids,t,p))[active]]
    for start in starts:
        try:
            sol=least_squares(residual,start,jac="3-point",xtol=1e-14,ftol=1e-14,gtol=1e-14,max_nfev=5000)
            assert sol.success,sol.message
            polish=root(residual,sol.x,method="hybr",options=dict(xtol=1e-12,maxfev=5000))
            if (polish.success or "xtol" in polish.message) and np.max(np.abs(residual(polish.x)))<np.max(np.abs(residual(sol.x))):
                sol.x=polish.x
            k=np.ones(len(z));k[active]=np.exp(sol.x)
            beta,x,y=a4.rachford_rice(z,k)
            fv,fl=a4.min_gibbs(model,t,p,y),a4.min_gibbs(model,t,p,x)
            eq=float(np.max(np.abs(np.log(y[active])+fv[2][active]-np.log(x[active])-fl[2][active])))
            assert eq<1e-11,("independent split fugacity",eq)
            assert 0<beta<1,("nonphysical split fraction",beta)
            assert fv[1]<fl[1]*.999,"trivial or nearly identical density solution"
            return beta,x,y,fv[0],fl[0],eq
        except (AssertionError,ValueError) as error:
            failures.append(str(error))
    raise AssertionError(("all independent split starts failed",failures))


def case(lib,case_id,ids,z,t,p):
    print("GERG reference",case_id,flush=True)
    z=np.array(z);model=Reference(lib,ids)
    model.trial_records=[]
    tm,w,zw=a4.tpd_search(model,ids,t,p,z)
    feed=a4.min_gibbs(model,t,p,z)
    base=dict(id=case_id,ids=ids,z=z.tolist(),T=t,p=p,feed_min_tm=tm)
    if tm>=-1e-8:
        props=model.properties(t,feed[1],z)
        return dict(**base,phases=1,composition=z.tolist(),**props,enthalpy=props["h"],reference_tpd_trajectories=model.trial_records)
    active=z>0
    ratio=np.zeros(len(z));ratio[active]=np.log(np.array(w)[active]/z[active])
    beta,x,y,zv,zl,eq=split(model,t,p,z,ratio if zw>feed[0] else -ratio)
    assert 0<beta<1
    phase_tm=[a4.tpd_search(model,ids,t,p,comp)[0] for comp in (x,y)]
    assert min(phase_tm)>=-1e-8,(case_id,phase_tm)
    vapor=a4.min_gibbs(model,t,p,y);liquid=a4.min_gibbs(model,t,p,x)
    assert vapor[1]<liquid[1]
    vp=model.properties(t,vapor[1],y);lp=model.properties(t,liquid[1],x)
    def g(w,r):
        return float(np.sum(w[active]*(np.log(w[active])+r[active])))
    dg=beta*g(y,vapor[2])+(1-beta)*g(x,liquid[2])-g(z,feed[2])
    assert dg<0
    return dict(**base,phases=2,vapor_fraction=float(beta),vapor=y.tolist(),liquid=x.tolist(),
        vapor_properties=vp,liquid_properties=lp,fugacity_residual=float(eq),phase_min_tm=phase_tm,
        reduced_gibbs_change=dg,enthalpy=beta*vp["h"]+(1-beta)*lp["h"],reference_tpd_trajectories=model.trial_records)


def saturation(lib,name,t):
    model=Reference(lib,[name]);x=np.array([1.])
    guesses=(12000.,250.) if name=="propane" else (15000.,1000.)
    liquid,vapor=model.model.pure_VLE_T(t,*guesses,50)
    assert liquid>vapor*1.01
    pressure=lambda d:d*R*t*(1+model.model.get_Ar01(t,d,x))
    p=pressure(vapor)
    assert abs(pressure(liquid)-p)<=1e-9*p
    assert abs(np.log(model.model.get_fugacity_coefficients(t,liquid*x))[0]-np.log(model.model.get_fugacity_coefficients(t,vapor*x))[0])<1e-11
    exterior=model.roots(t,p,x)
    assert abs(exterior[0][1]-vapor)<2e-5 and abs(exterior[-1][1]-liquid)<2e-5
    return dict(component=name,T=t,p=p,vapor_properties=model.properties(t,vapor,x),liquid_properties=model.properties(t,liquid,x))


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output",type=Path,required=True)
    parser.add_argument("--case",help="Generate one diagnostic candidate case only")
    args=parser.parse_args()
    assert teqp.__version__=="0.23.2"
    with tempfile.TemporaryDirectory(prefix="gerg-equilibrium-reference-") as folder:
        tmp=Path(folder);(tmp/"reference.cpp").write_text(CODE)
        subprocess.run(["g++","-O2","-std=c++11","-ftrivial-auto-var-init=zero","-shared","-fPIC","-I",str(SOURCE),str(tmp/"reference.cpp"),"-o",str(tmp/"reference.so")],check=True)
        lib=ctypes.CDLL(str(tmp/"reference.so"))
        lib.setup()
        lib.density_roots.argtypes=[ctypes.c_double,ctypes.c_double,PTR,PTR];lib.density_roots.restype=ctypes.c_int
        lib.properties.argtypes=[ctypes.c_double,ctypes.c_double,PTR,PTR];lib.properties.restype=None
        a4.roots=lambda model,t,p,x:model.roots(t,p,x)
        a4.tpd_search=tpd_search
        c2_sat=saturation(lib,"ethane",260.)
        inputs=[
            ("binary-split",["methane","propane"],[.5,.5],250.,3e6),
            ("binary-vapor",["methane","propane"],[.5,.5],330.,2e6),
            ("binary-liquid",["methane","propane"],[.05,.95],250.,5e6),
            ("binary-zero-n2",["methane","propane","nitrogen"],[.5,.5,0.],250.,3e6),
            ("ng-vapor",a4.NG,[.88,.06,.03,.02,.01],330.,5e6),
            ("ng-rich-split",a4.NG,[.70,.15,.13,.01,.01],250.,3e6),
            ("ng-rich-low-T",a4.NG,[.70,.15,.13,.01,.01],230.,2e6),
            ("ng-mixed",a4.NG,a4.MIX,240.,5e6),
            ("binary-boundary-low",["methane","propane"],[.5,.5],290.,1e6),
            ("binary-boundary-high",["methane","propane"],[.5,.5],290.,5e6),
            ("pure-propane",["propane"],[1.],300.,1e6),
            ("pure-methane",["methane"],[1.],300.,5e6),
            ("narrow-ethane-propane",["ethane","propane"],[.98,.02],260.,.97*c2_sat["p"]),
        ]
        if args.case:
            inputs=[c for c in inputs if c[0]==args.case]
            assert inputs,"unknown diagnostic case"
        cases=[case(lib,*c) for c in inputs]
        ph=[]
        if not args.case:
            for c in cases:
                if c["id"] in ["binary-split","binary-vapor","binary-liquid","ng-rich-split","ng-vapor","narrow-ethane-propane"]:
                    ph.append(dict(id=c["id"]+"-ph",ids=c["ids"],z=c["z"],p=c["p"],target_h=c["enthalpy"],T=c["T"],
                        temperature_bracket=[max(230.,c["T"]-15),min(350.,c["T"]+15)],expected_tp=c))
            for t in [280.,300.]:
                sat=saturation(lib,"propane",t)
                for beta in [0.,.25,.75,1.]:
                    h=beta*sat["vapor_properties"]["h"]+(1-beta)*sat["liquid_properties"]["h"]
                    ph.append(dict(id=f"pure-propane-{t}-{beta}",ids=["propane"],z=[1.],p=sat["p"],target_h=h,T=t,
                        temperature_bracket=[250.,330.],vapor_fraction=beta,saturation=sat))
            inlet=next(c for c in cases if c["id"]=="ng-vapor")
            ids,z,p=a4.NG,np.array(inlet["z"]),3e6
            valve=Reference(lib,ids);valve.trial_records=[]
            def energy(t):
                r=a4.min_gibbs(valve,t,p,z)
                return valve.properties(t,r[1],z)["h"]-inlet["enthalpy"]
            t=brentq(energy,300.,350.,xtol=1e-10,rtol=2e-14)
            outlet=case(lib,"ng-isenthalpic-valve-outlet",ids,z.tolist(),t,p)
            assert abs(outlet["enthalpy"]-inlet["enthalpy"])<=2e-6
            ph.append(dict(id="ng-isenthalpic-valve",ids=ids,z=z.tolist(),p=p,target_h=inlet["enthalpy"],T=t,
                temperature_bracket=[300.,350.],expected_tp=outlet,inlet= dict(T=inlet["T"],p=inlet["p"],h=inlet["enthalpy"])))
    output=dict(schema="gerg-equilibrium-v1",packages={n:importlib.metadata.version(n) for n in ["teqp","numpy","scipy"]},
        source_sha256=hashlib.sha256((SOURCE/"GERG2008.cpp").read_bytes()).hexdigest(),
        a4_adapter_sha256=hashlib.sha256((ROOT/"validation/generate-flash-references.py").read_bytes()).hexdigest(),
        generator_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),seed=a4.SEED,random_tpd_starts=a4.RANDOM_STARTS,
        independence="NIST native pressure/calorics, teqp AD fugacity, SciPy L-BFGS-B TPD and hybrid Powell split; no Rust calls; outer roots from denser independent sampling",
        scope="dry five-component fluid EOS, regular finite-start VLE evidence; no solids, global stability proof, critical-region or experimental accuracy claim",
        tolerances=dict(composition=[2e-8,1e-7],vapor_fraction=[2e-8,0.],temperature=[2e-6,0.],rho=[2e-5,2e-8],Z=[2e-9,2e-8],h=[2e-5,2e-9],s=[2e-7,2e-9],cp=[2e-5,2e-8]),tp=cases,ph=ph)
    args.output.write_text(json.dumps(output,indent=2,allow_nan=False)+"\n")
    print(f"Wrote {len(cases)} GERG TP and {len(ph)} PH cases to {args.output}")


if __name__=="__main__":main()
