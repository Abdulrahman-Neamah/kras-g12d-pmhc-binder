"""Does the published high-discrimination TCR contact the G12D Asp side chain?

7OW6 = JDIa41b1 TCR on HLA-A*11:01 / VVVGADGVGK (G12D)
7OW5 = the same TCR on VVVGAGGVGK (wild type)
chains: A = HLA heavy, B = b2m, C = peptide, D = TCRalpha, E = TCRbeta
"""
import sys, warnings; warnings.filterwarnings('ignore')
import numpy as np
from Bio.PDB import PDBParser, ShrakeRupley
from scipy.spatial import distance

import os
# Directory holding 7OW5.pdb and 7OW6.pdb (download from the RCSB PDB).
# Override with:  STRUCT_DIR=/path/to/pdbs python3 01_position6_readout.py
D = os.environ.get('STRUCT_DIR', './structures') + '/'
BB = {'N', 'CA', 'C', 'O', 'OXT'}


def load(pid):
    return PDBParser(QUIET=True).get_structure(pid, D + pid + '.pdb')[0]


def atoms(ch, want_het=False):
    return [a for r in ch for a in r
            if (r.id[0] == ' ') or (want_het and r.get_resname() == 'HOH')]


def report(pid):
    m = load(pid)
    pep = m['C']
    r6 = pep[6]
    print('=' * 78)
    print(f'{pid}   peptide position 6 = {r6.get_resname()}')
    sc = [a for a in r6 if a.get_id() not in BB]
    bb = [a for a in r6 if a.get_id() in BB]
    print(f'  side-chain atoms modelled: {[a.get_id() for a in sc] or "none (glycine)"}')

    tcr = atoms(m['D']) + atoms(m['E'])
    tcr_lab = [f'{a.get_parent().get_parent().id}/{a.get_parent().get_resname()}'
               f'{a.get_parent().id[1]} {a.get_id()}' for a in tcr]
    tc = np.array([a.coord for a in tcr])

    # --- side chain vs TCR
    if sc:
        Dm = distance.cdist(np.array([a.coord for a in sc]), tc)
        i, j = np.unravel_index(np.argmin(Dm), Dm.shape)
        print(f'  SIDE CHAIN -> nearest TCR atom: {Dm.min():.2f} A '
              f'({sc[i].get_id()} -- {tcr_lab[j]})')
        for cut in (4.0, 5.0, 6.0):
            print(f'     TCR atoms within {cut} A of the side chain: {(Dm < cut).sum()}')

    # --- backbone vs TCR
    Db = distance.cdist(np.array([a.coord for a in bb]), tc)
    i, j = np.unravel_index(np.argmin(Db), Db.shape)
    print(f'  BACKBONE  -> nearest TCR atom: {Db.min():.2f} A '
          f'({bb[i].get_id()} -- {tcr_lab[j]})')
    for cut in (4.0, 5.0):
        print(f'     TCR atoms within {cut} A of the backbone: {(Db < cut).sum()}')

    # --- side chain vs HLA
    if sc:
        hla = atoms(m['A'])
        hl = np.array([a.coord for a in hla])
        Dh = distance.cdist(np.array([a.coord for a in sc]), hl)
        print('  SIDE CHAIN -> HLA polar partners under 3.6 A:')
        for ii, jj in np.argwhere(Dh < 3.6):
            a, b = sc[ii], hla[jj]
            if a.element in 'NO' and b.element in 'NO':
                print(f'     {a.get_id():4s} -- A/{b.get_parent().get_resname()}'
                      f'{b.get_parent().id[1]} {b.get_id():4s}  {Dh[ii, jj]:.2f} A')

    # --- burial of position 6 by HLA vs by TCR
    sr = ShrakeRupley()
    import copy
    def sasa_of(chains):
        sub = copy.deepcopy(m)
        for ch in list(sub):
            if ch.id not in chains:
                sub.detach_child(ch.id)
        for ch in sub:
            for r in list(ch):
                if r.id[0] != ' ':
                    ch.detach_child(r.id)
        sr.compute(sub, level='R')
        r = sub['C'][6]
        tot = r.sasa
        side = sum(a.sasa for a in r if a.get_id() not in BB)
        return tot, side
    free_t, free_s = sasa_of({'C'})
    mhc_t, mhc_s = sasa_of({'A', 'B', 'C'})
    all_t, all_s = sasa_of({'A', 'B', 'C', 'D', 'E'})
    print(f'  SASA of residue 6 (total / side chain), A^2:')
    print(f'     peptide alone      {free_t:7.1f} {free_s:7.1f}')
    print(f'     + HLA (pMHC)       {mhc_t:7.1f} {mhc_s:7.1f}   buried by HLA: '
          f'{free_t - mhc_t:.1f} / {free_s - mhc_s:.1f}')
    print(f'     + HLA + TCR        {all_t:7.1f} {all_s:7.1f}   buried by TCR: '
          f'{mhc_t - all_t:.1f} / {mhc_s - all_s:.1f}')
    return m


m6 = report('7OW6')
m5 = report('7OW5')
