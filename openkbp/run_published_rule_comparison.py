#!/usr/bin/env python3
"""C4 published-rule comparison for the 60 Experiment-A OpenKBP controlled lesions.

Published methods:
- Ferreira et al. 2015 method TV: TV = 95% prescription isodose; InField >95% lesion in TV; Marginal 20-95%; OutOfField <20%.
- Mohamed et al. 2016 combined centroid/fD95 typology: native A-E labels retained.

Li et al. 2014 is not approximated because OpenKBP lacks its required GTV/CTV1/CTV2 hierarchy.
"""
from pathlib import Path
import argparse,re
import numpy as np
import pandas as pd
from scipy.stats import beta
SHAPE=(128,128,128); NVOX=np.prod(SHAPE)
KNOWN=["PTV70","PTV63","PTV56","Brainstem","SpinalCord","RightParotid","LeftParotid","Esophagus","Larynx","Mandible"]
def read_sparse_csv(path,binary=False):
    arr=np.zeros(NVOX,dtype=np.float32); df=pd.read_csv(path,index_col=0)
    raw=pd.to_numeric(pd.Series(df.index),errors="coerce"); valid=raw.notna().values; idx=raw[valid].astype(np.int64).values
    vals=np.ones(len(idx),dtype=np.float32) if binary or "data" not in df.columns else pd.to_numeric(df["data"],errors="coerce").fillna(0).values.astype(np.float32)[valid]
    good=(idx>=0)&(idx<NVOX); arr[idx[good]]=vals[good]; return arr.reshape(SHAPE)
def spacing(path):
    x=pd.to_numeric(pd.Series(pd.read_csv(path,header=None).values.ravel()),errors="coerce").dropna().values.astype(float); return x[-3:]
def load(p):
    s={}; p=Path(p)
    for name in KNOWN:
        f=p/f"{name}.csv"
        if f.exists(): s[name]=read_sparse_csv(f,True).astype(bool)
    return {"id":p.name,"dose":read_sparse_csv(p/"dose.csv"),"possible":read_sparse_csv(p/"possible_dose_mask.csv",True).astype(bool),"spacing":spacing(p/"voxel_dimensions.csv"),"structures":s}
def rx(s):
    v=[float(m.group(1)) for x in s if (m:=re.fullmatch(r"PTV(\d+)",x))]; return max(v) if v else 70.
def sphere(shape,sp,c,r=10.):
    sp=np.asarray(sp,float); c=np.asarray(c,float); rv=np.ceil(r/sp).astype(int); lo=np.maximum(np.floor(c-rv-1).astype(int),0); hi=np.minimum(np.ceil(c+rv+1).astype(int),np.asarray(shape)-1)
    I,J,K=np.meshgrid(np.arange(lo[0],hi[0]+1),np.arange(lo[1],hi[1]+1),np.arange(lo[2],hi[2]+1),indexing="ij")
    d2=((I-c[0])*sp[0])**2+((J-c[1])*sp[1])**2+((K-c[2])*sp[2])**2; m=np.zeros(shape,bool); m[lo[0]:hi[0]+1,lo[1]:hi[1]+1,lo[2]:hi[2]+1]=d2<=r*r; return m
def seed_center(dose,R,possible,s,kind):
    n=dose/R; ptv=s.get("PTV70",s.get("PTV63",s.get("PTV56",possible)))
    if kind=="high-dose-seeded": cand=possible&ptv&(n>=.95); target=1.
    elif kind=="gradient-seeded": cand=possible&(n>=.35)&(n<.90); target=.65
    else: cand=possible&(n<.20); target=.10
    idx=np.argwhere(cand)
    if len(idx)==0: e=np.abs(n-target); e[~possible]=np.inf; return np.array(np.unravel_index(np.argmin(e),dose.shape))
    ctr=np.asarray(dose.shape)/2; return idx[np.argmin(np.sum((idx-ctr)**2,axis=1))]
def metrics(dose,m,R):
    v=dose[m]; q=v/R; return {"f_high":float(np.mean(q>=.95)),"f_mid":float(np.mean((q>=.20)&(q<.95))),"f_low":float(np.mean(q<.20)),"d95":float(np.percentile(v,5)),"dmean":float(v.mean())}
def dominant(m): return max({"InField":m["f_high"],"Marginal":m["f_mid"],"OutOfField":m["f_low"]},key=lambda x:{"InField":m["f_high"],"Marginal":m["f_mid"],"OutOfField":m["f_low"]}[x])
def strict(m): return "InField" if m["f_high"]>=.95 else ("OutOfField" if m["f_low"]>=.95 else "Marginal")
def ferreira(m): return "InField" if m["f_high"]>.95 else ("Marginal" if m["f_high"]>=.20 else "OutOfField")
def origin(c,s):
    ijk=tuple(int(x) for x in c)
    for name,R in [("PTV70",70.),("PTV63",63.),("PTV56",56.)]:
        if name in s and s[name][ijk]: return name,R
    return None,None
def mohamed(c,m,s):
    o,R=origin(c,s)
    if o is None:return "E",o,R
    covered=m["d95"]>=.95*R
    return (("A" if covered else "B") if o=="PTV70" else ("C" if covered else "D")),o,R
def cp(k,n): return (0. if k==0 else float(beta.ppf(.025,k,n-k+1)),1. if k==n else float(beta.ppf(.975,k+1,n-k)))
def main(root,out):
    root=Path(root); out=Path(out); out.mkdir(parents=True,exist_ok=True); pts=sorted([p for p in root.rglob("pt_*") if (p/"dose.csv").exists()],key=lambda p:int(re.findall(r"\d+",p.name)[-1])); rows=[]
    for pp in pts:
        p=load(pp); R=rx(p["structures"])
        for kind in ["high-dose-seeded","gradient-seeded","low-dose-seeded"]:
            c=seed_center(p["dose"],R,p["possible"],p["structures"],kind); m=metrics(p["dose"],sphere(p["dose"].shape,p["spacing"],c),R); mt,o,oR=mohamed(c,m,p["structures"])
            d=dominant(m); st=strict(m); f=ferreira(m); rows.append({"patient":p["id"],"seed_type":kind,**m,"dominant_fraction_v1":d,"strict_overlap_v1":st,"ferreira_tv":f,"mohamed_type":mt,"mohamed_origin_ptv":o,"dominant_vs_ferreira_disagree":d!=f,"strict_vs_ferreira_disagree":st!=f})
    df=pd.DataFrame(rows); df.to_csv(out/"experiment_A_published_rule_comparison.csv",index=False)
    expected={"gradient-seeded":0,"high-dose-seeded":19,"low-dose-seeded":2}; observed={s:int(((g:=df[df.seed_type==s]).dominant_fraction_v1!=g.strict_overlap_v1).sum()) for s in expected}
    if observed!=expected: raise RuntimeError(f"Experiment-A reconstruction mismatch: {observed}")
    summary=[]
    for seed,g in df.groupby("seed_type"):
        for name,col in [("dominant_fraction_v1 vs Ferreira-TV","dominant_vs_ferreira_disagree"),("strict_overlap_v1 vs Ferreira-TV","strict_vs_ferreira_disagree")]:
            k=int(g[col].sum()); lo,hi=cp(k,len(g)); summary.append([seed,name,k,len(g),k/len(g),lo,hi])
    pd.DataFrame(summary,columns=["seed_type","comparison","disagreements","n","fraction","ci95_low","ci95_high"]).to_csv(out/"published_rule_pairwise_summary.csv",index=False)
    print("reconstruction",observed); print(pd.crosstab(df.seed_type,df.ferreira_tv)); print(pd.crosstab(df.seed_type,df.mohamed_type))
if __name__=="__main__":
    a=argparse.ArgumentParser(); a.add_argument("openkbp_root"); a.add_argument("--out",default="c4_results"); x=a.parse_args(); main(x.openkbp_root,x.out)
