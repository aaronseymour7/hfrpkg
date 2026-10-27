"""Recompute every benchmark error (data/benchmark_errors.csv) from the exported energies in data/si/
and the reactions in data/objective_check.csv, and list the ones that differ by more than 0.01 kcal/mol.
With --write, data/benchmark_errors.csv is rewritten with the recomputed errors (used after the
saddle-point species were reoptimized, and to use the lowest-energy copy of species computed in two runs).

    cd paper && python recompute_from_structures.py [--write]
"""
import csv,glob,re,sys
from rdkit import Chem, RDLogger
RDLogger.DisableLog('rdApp.*')
H2K=2625.4996/4.184
LV=[('b3lyp','B3LYP/6-31+G(d)'),('m06','M06-2X/6-31+G(d)'),('b97','wB97X-D/def2-TZVP')]
# every species, in every reaction, uses its lowest-energy copy (E + ZPVE) across the exports; on a tie the
# later run (structures_2, the new-species runs) is used
ES={}
for f in sorted(glob.glob('data/si/*.xyz')):
    L=open(f).read().splitlines();i=0
    while i<len(L) and L[i].strip():
        k=int(L[i]);p=L[i+1].split(' ');kv=dict(re.findall(r'(\S+?)=(\S+)',' '.join(p[4:])))
        lv=[n for r,n in LV if r in p[0].lower()][0]
        e=(float(kv['E']),float(kv['ZPVE']),float(kv['E_sp']) if 'E_sp' in kv else None)
        if (lv,p[3]) not in ES or e[0]+e[1]<=sum(ES[(lv,p[3])][:2]): ES[(lv,p[3])]=e
        i+=2+k
atct={}
for l in open('../hfrpkg/data/ATcT_lib.txt'):
    t=l.rstrip('\n').split('\t')
    if len(t)>=6: atct[t[3]]=float(t[5])
inch=lambda s: Chem.MolToInchi(Chem.MolFromSmiles(s))
def side(x):
    out=[]
    for t in x.split(' + '):
        t=t.strip().split(' ')
        out.append((int(t[0]),t[1]) if len(t)==2 else (1,t[0]))
    return out
rx={}
for r in csv.DictReader(open('data/objective_check.csv')):
    rx[(r['smiles'],r['scheme'])]=r['legacy_no_imine'] if r['scheme']=='HY' else r['current']
# reactions actually run in the benchmark where they differ from the regenerated ones (index.txt of the run folders)
rx[('NCC(=O)O','HO')]='CC + CO + NCC(=O)O -> CN + CC(=O)O + CCO'
A, B = 'OCl + OOOCl -> 2 OOCl', 'OO + OOOCl -> OOCl + OOO'  # HClO3: the two runs differ by level
LEVEL_RX = {('OOOCl', 'HY', 'B3LYP/6-31+G(d)'): A, ('OOOCl', 'HO', 'B3LYP/6-31+G(d)'): A,
            ('OOOCl', 'HY', 'M06-2X/6-31+G(d)'): A, ('OOOCl', 'HO', 'M06-2X/6-31+G(d)'): B,
            ('OOOCl', 'HY', 'wB97X-D/def2-TZVP'): B, ('OOOCl', 'HO', 'wB97X-D/def2-TZVP'): A}
rx[('O=c1cc1','HO')]='CC=O + C=O + C=CC + O=c1cc1 -> 2 C=CC=O + CC(C)=O'
worst=0;n=0;bad=[]
rows=list(csv.DictReader(open('data/benchmark_errors.csv')))
for r in rows:
    for sch in ('ID','HY','HO'):
        lhs,rhs=LEVEL_RX.get((r['smiles'],sch,r['geometry_level']),rx[(r['smiles'],sch)]).split(' -> ')
        tot=0;hf=0
        for sign,sd in ((-1,side(lhs)),(1,side(rhs))):
            for nu,s in sd:
                ii=inch(s);e,z,sp=ES[(r['geometry_level'],ii)]
                tot+=sign*nu*((sp if r['stage']=='sp' else e)+z)
                if ii!=r['inchi']: hf+=sign*nu*atct[ii]
        # target appears on lhs with coefficient 1
        calc=hf/4.184 - tot*H2K
        err=calc-float(r['atct_kcal']); d=abs(err-float(r[sch])); n+=1
        if d>0.01: bad.append((r['smiles'],r['stage'],r['geometry_level'],sch,round(err,3),round(float(r[sch]),3)))
        worst=max(worst,d)
        r['new_'+sch]=err
print(n,'errors recomputed; max diff',round(worst,4),'; >0.01:',len(bad))
for b in bad: print(b)
if '--write' in sys.argv:
    cols=list(rows[0].keys())[:9]
    with open('data/benchmark_errors.csv','w',newline='') as fh:
        w=csv.writer(fh); w.writerow(cols)
        for r in rows: w.writerow([r[c] if c not in ('ID','HY','HO') else repr(r['new_'+c]) for c in cols])
    print('data/benchmark_errors.csv rewritten')
