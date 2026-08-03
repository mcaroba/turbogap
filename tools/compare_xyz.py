#!/usr/bin/env python3
"""Compare energies/forces/virial between two TurboGAP extended-XYZ trajectories."""
import sys, re
import numpy as np

def parse_frames(path):
    frames=[]
    with open(path) as f:
        while True:
            line=f.readline()
            if not line: break
            if not line.strip(): continue
            n=int(line.split()[0])
            comment=f.readline()
            props=dict(re.findall(r'(\w+)=("[^"]*"|\S+)',comment))
            cols=[]
            for _ in range(n):
                cols.append(f.readline().split())
            frames.append((props,cols))
    return frames

def getf(props,key):
    v=props.get(key)
    return None if v is None else float(v.strip('"'))

def getarr(props,key):
    v=props.get(key)
    return None if v is None else np.array([float(x) for x in v.strip('"').split()])

def main(a,b,ftol=1e-4,etol=1e-6):
    A=parse_frames(a); B=parse_frames(b)
    print(f"frames: {a}={len(A)}  {b}={len(B)}")
    nf=min(len(A),len(B))
    worst=0.0; fail=False
    for i in range(nf):
        pa,ca=A[i]; pb,cb=B[i]
        print(f"\n--- frame {i} ---")
        for k in ["energy","energy_soap","energy_2b","energy_3b","energy_core_pot","energy_vdw"]:
            va,vb=getf(pa,k),getf(pb,k)
            if va is None or vb is None: continue
            d=abs(va-vb); rel=d/max(abs(va),1e-30)
            flag="OK " if rel<1e-6 or d<1e-6 else "DIFF"
            print(f"  {flag} {k:16s} {va:20.8f} {vb:20.8f}  absdiff={d:.3e} rel={rel:.2e}")
            if flag=="DIFF": fail=True
        for k in ["virial"]:
            va,vb=getarr(pa,k),getarr(pb,k)
            if va is None or vb is None: continue
            d=np.abs(va-vb).max(); rel=d/max(np.abs(va).max(),1e-30)
            flag="OK " if rel<1e-6 else "DIFF"
            print(f"  {flag} {k:16s} maxabsdiff={d:.6e} rel={rel:.2e}")
            print(f"       cpu={va}")
            print(f"       gpu={vb}")
            if flag=="DIFF": fail=True
        # forces: species,x,y,z,fx,fy,fz,...
        fa=np.array([[float(x) for x in r[4:7]] for r in ca])
        fb=np.array([[float(x) for x in r[4:7]] for r in cb])
        if fa.shape==fb.shape:
            d=np.abs(fa-fb); md=d.max()
            fmag=np.abs(fa).max()
            worst=max(worst,md)
            idx=np.unravel_index(d.argmax(),d.shape)
            flag="OK " if md<ftol else "DIFF"
            print(f"  {flag} forces           maxabsdiff={md:.6e} (|F|max={fmag:.4f}) rms={np.sqrt((d**2).mean()):.3e}")
            print(f"       worst atom {idx[0]} comp {idx[1]}: cpu={fa[idx]:.8f} gpu={fb[idx]:.8f}")
            if flag=="DIFF": fail=True
        # local energies col 7
        try:
            la=np.array([float(r[7]) for r in ca]); lb=np.array([float(r[7]) for r in cb])
            d=np.abs(la-lb).max()
            print(f"  {'OK ' if d<1e-6 else 'DIFF'} local_energy     maxabsdiff={d:.6e}")
        except Exception: pass
    print("\nRESULT:", "FAIL" if fail else "PASS")
    return 1 if fail else 0

if __name__=="__main__":
    sys.exit(main(sys.argv[1],sys.argv[2]))
