#!/usr/bin/env python3
"""Assemble the Challenge 5 submission CSV from ipsae_scores.csv.

  python3 build_submission.py ipsae_scores.csv out.csv [--max 20]

Ranks by the worse of the two discrimination gaps (WT and apo), because a design must
clear both goals, not average them. Drops anything with no real G12D interface.
"""
import csv, sys, statistics as st

MIN_G12D   = 0.15   # below this there is no interface to speak of
SEQS = {}
# Sequence sources. Design sequences are withheld from this repository pending
# competition results, so these files are not present here; point SEQ_DIR at your own.
import os
_d = os.environ.get('SEQ_DIR', '.')
for p in (f'{_d}/ch5_mpnn_variants.csv', f'{_d}/ch5_d010_extra.csv', f'{_d}/ch5_candidates.csv'):
    try:
        for r in csv.DictReader(open(p)):
            n = r.get('name') or r['design_id'].replace('rfd4_design_', 'rfd4_d')
            SEQS.setdefault(n, r['sequence'])
    except FileNotFoundError:
        pass

src, out = sys.argv[1], sys.argv[2]
cap = int(sys.argv[sys.argv.index('--max') + 1]) if '--max' in sys.argv else 20
f = lambda v: float('nan') if v in ('', 'nan', None) else float(v)
base = lambda n: n.rsplit('_s', 1)[0] if '_s' in n and n.rsplit('_s', 1)[1].isdigit() else n

G = {}
for r in csv.DictReader(open(src)):
    G.setdefault(base(r['name']), {}).setdefault(r['target'], []).append(r)
m = lambda L: st.mean(L) if L else float('nan')

rows = []
for d, t in G.items():
    gp = [f(x['ipsae_peptide_max']) for x in t.get('G12D', [])]
    wp = [f(x['ipsae_peptide_max']) for x in t.get('WT', [])]
    gh = [f(x['ipsae_HLA_max'])     for x in t.get('G12D', [])]
    ah = [f(x['ipsae_HLA_max'])     for x in t.get('apo', [])]
    if not gp or m(gp) < MIN_G12D:
        continue
    wgap = m(gp) - m(wp) if wp else float('nan')
    agap = m(gh) - m(ah) if ah else float('nan')
    both = [x for x in (wgap, agap) if x == x]
    rows.append({'name': d, 'seeds': len(gp), 'g12d': m(gp), 'wt_gap': wgap,
                 'apo_gap': agap, 'worst': min(both) if both else float('-inf'),
                 'seq': SEQS.get(d, '')})

rows.sort(key=lambda r: -r['worst'])
missing = [r['name'] for r in rows if not r['seq']]
if missing:
    print('NO SEQUENCE FOUND for:', missing)

keep = rows[:cap]
with open(out, 'w', newline='') as fh:
    w = csv.writer(fh)
    w.writerow(['name', 'sequence', 'molecule_class'])
    for r in keep:
        if r['seq']:
            w.writerow([r['name'], r['seq'], 'protein'])

print(f"{'name':18s} {'n':>1s} {'G12D':>6s} {'WTgap':>7s} {'apoGap':>7s} {'worst':>7s} {'len':>4s}")
for r in keep:
    s = lambda x: '   -  ' if x != x else f'{x:6.3f}'
    print(f"{r['name']:18s} {r['seeds']:1d} {s(r['g12d']):>6s} {s(r['wt_gap']):>7s} "
          f"{s(r['apo_gap']):>7s} {s(r['worst']):>7s} {len(r['seq']):4d}")
print(f'\nwrote {len([r for r in keep if r["seq"]])} designs -> {out}')
