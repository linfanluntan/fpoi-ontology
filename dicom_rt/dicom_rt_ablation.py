#!/usr/bin/env python3
from pathlib import Path
import argparse,shutil,pydicom
def modality(p):
    try:return str(pydicom.dcmread(str(p),stop_before_pixels=True,force=True).Modality)
    except:return ""
def copy_without(src,dst,rm):
    if dst.exists():shutil.rmtree(dst)
    dst.mkdir(parents=True); n=0
    for f in src.rglob("*"):
        if not f.is_file():continue
        if modality(f)==rm:n+=1;continue
        o=dst/f.relative_to(src);o.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(f,o)
    return n
if __name__=="__main__":
    ap=argparse.ArgumentParser();ap.add_argument("package");ap.add_argument("--out",default="dicom_rt_ablations");a=ap.parse_args()
    for m in ["CT","RTSTRUCT","RTPLAN","RTDOSE"]:print(m,copy_without(Path(a.package),Path(a.out)/("minus_"+m),m))
