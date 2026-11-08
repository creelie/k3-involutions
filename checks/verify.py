"""Exact finite combinatorics supporting the geometric proofs; no numerical geometry."""
import json
from collections import Counter
from pathlib import Path

ROOT=Path(__file__).resolve().parent
def pair(c): return (c >> 2, c & 3)
def dot(c,d): return (c & d).bit_count() % 2
def trace(c,r):
    # Independent source: node permutation plus Kunneth action on H2(A).
    node_trace=sum((v ^ c)==v for v in range(16))
    return node_trace + (6 if not r else -2)
def kind(c,r):
    a,b=pair(c)
    if not r: return 'identity' if c==0 else 'eight_points'
    if a and b: return 'free'
    if a: return 'two_F'
    if b: return 'two_E'
    return 'eight_rational'

elements=[(c,r) for c in range(16) for r in range(2)]
assert len(elements)==32
for c,r in elements:
    assert (c^c,r^r)==(0,0)
    for d,s in elements: assert (c^d,r^s) in elements
counts=Counter(kind(c,r) for c,r in elements if (c,r)!=(0,0))
assert counts=={'eight_points':15,'free':9,'two_E':3,'two_F':3,'eight_rational':1}
multiplicities={}
for d in range(16):
    for s in range(2):
        total=sum(trace(c,r)*(-1)**(dot(c,d)+s*r) for c,r in elements)
        assert total%32==0
        mult=total//32
        expected=(3 if s==0 else 4) if d==0 else (1 if s==0 else 0)
        assert mult==expected
        multiplicities[f'{d:04b}:{s}']=mult
assert sum(multiplicities.values())==22

descents=[]
for c0 in range(16):
    a,b=pair(c0)
    if not(a and b): continue
    dc=Counter()
    for c in range(1,16):
        x,y=pair(c)
        if c==c0: k='four_rational'
        elif (x,y)==(a,0): k='two_E_quotients'
        elif (x,y)==(0,b): k='two_F_quotients'
        elif x==a and y not in (0,b): k='one_E'
        elif y==b and x not in (0,a): k='one_F'
        else: k='no_curve'
        dc[k]+=1
        tr=(trace(c,0)+trace(c^c0,1))//2
        assert tr==(10 if c==c0 else 2)
        euler=12 if k=='four_rational' else 4
        assert euler==2+tr
        # Character restriction gives the same trace without the Lefschetz formula.
        direct=2+sum((-1)**dot(c,d) for d in range(16) if dot(c0,d)==0)
        assert direct==tr
    assert dc=={'four_rational':1,'two_E_quotients':1,'two_F_quotients':1,'one_E':2,'one_F':2,'no_curve':8}
    descents.append({'c0':pair(c0),'types':dict(dc),'numerical_kernel':[0,c0]})

# Exact gcd certificate for the Fermat root polynomial and its derivative:
# (z^4+1) - (z/4)(4z^3) = 1.
f=[1,0,0,0,1]
product=[0,0,0,0,1]
assert [a-b for a,b in zip(f,product)]==[1,0,0,0,0]
result={'status':'passed','arithmetic':'integers and finite bit operations',
 'group_order':32,'nonidentity_fixed_locus_counts':dict(counts),
 'H2_character_multiplicities':multiplicities,'Enriques_descents':descents,
 'Fermat_fixed_points':8,
 'scope':'Checks finite labels and characters. The geometric proof is in main.tex. No Lean formalization or novelty certification.'}
(ROOT/'results.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
print(json.dumps({'status':'passed','group_order':32,'descents':len(descents),'Fermat_fixed_points':8}))
