"""Which peptide residues does the validated TCR actually engage, and how?"""
import warnings; warnings.filterwarnings('ignore')
import numpy as np, copy
from Bio.PDB import PDBParser, ShrakeRupley
from scipy.spatial import distance
import os
# Directory holding 7OW5.pdb and 7OW6.pdb (download from the RCSB PDB).
# Override with:  STRUCT_DIR=/path/to/pdbs python3 02_tcr_footprint.py
D = os.environ.get('STRUCT_DIR', './structures') + '/'
BB={'N','CA','C','O','OXT'}
sr=ShrakeRupley()

def run(pid):
    m=PDBParser(QUIET=True).get_structure(pid,D+pid+'.pdb')[0]
    tcr=[a for c in ('D','E') for r in m[c] if r.id[0]==' ' for a in r]
    tc=np.array([a.coord for a in tcr])
    def sasa(chains):
        sub=copy.deepcopy(m)
        for ch in list(sub):
            if ch.id not in chains: sub.detach_child(ch.id)
        for ch in sub:
            for r in list(ch):
                if r.id[0]!=' ': ch.detach_child(r.id)
        sr.compute(sub, level='A')
        return {r.id[1]: (sum(a.sasa for a in r), sum(a.sasa for a in r if a.get_id() not in BB)) for r in sub['C']}
    pm=sasa({'A','B','C'}); full=sasa({'A','B','C','D','E'})
    print('='*86)
    print(f'{pid}  peptide {"".join(r.get_resname()[0] for r in m["C"])}   TCR footprint on the peptide')
    print(f'{"pos":>3} {"res":>4} {"SASA in pMHC":>13} {"buried by TCR":>14} {"sc buried":>10} {"min dist":>9}  nearest TCR atom')
    tot_bsa=0
    for r in m['C']:
        n=r.id[1]
        dm=distance.cdist(np.array([a.coord for a in r]), tc)
        i,j=np.unravel_index(np.argmin(dm),dm.shape)
        bsa=pm[n][0]-full[n][0]; scb=pm[n][1]-full[n][1]; tot_bsa+=bsa
        tag=f'{tcr[j].get_parent().get_parent().id}/{tcr[j].get_parent().get_resname()}{tcr[j].get_parent().id[1]} {tcr[j].get_id()}'
        star=' <<<' if bsa>15 else ''
        print(f'{n:3d} {r.get_resname():>4} {pm[n][0]:10.1f} A2 {bsa:11.1f} A2 {scb:8.1f} {dm.min():8.2f} A  {tag}{star}')
    print(f'    total peptide surface buried by TCR: {tot_bsa:.0f} A2')
run('7OW6'); run('7OW5')
