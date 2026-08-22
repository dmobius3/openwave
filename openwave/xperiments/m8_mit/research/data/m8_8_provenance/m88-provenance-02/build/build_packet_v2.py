import sys; sys.path.insert(0,'.')
from complex import *
import json, hashlib, pathlib
d=json.load(open('m3_v2.json')); s_id,t_id=d['s_id'],d['t_id']
M3=[{int(k):v for k,v in c.items()} for c in d['M3']]
R1=[(0,1)]*3+[(1,-1),(0,-1)]*2                       # s^3 (st)^-2
R2=[(1,1)]*5+[(1,-1),(0,-1)]*2                       # t^5 (st)^-2
M2=[[fox(r,k,[s_id,t_id]) for k in (0,1)] for r in (R1,R2)]
M1=[[rz((1,s_id),(-1,E))],[rz((1,t_id),(-1,E))]]
enc=lambda x: [[c,g] for g,c in sorted(x.items()) if c]
GP=pathlib.Path('m8_5a_packet.json')
src=b"".join(pathlib.Path(f).read_bytes() for f in
             ["qphi.py","complex.py","kernel.py","sat.py","search_sym.py","build_packet_v2.py"])
packet={
 "format_version":"m8_8-construction-1",
 "group_packet_sha256":hashlib.sha256(GP.read_bytes()).hexdigest(),
 "abstract_generators":{"s":s_id,"t":t_id},
 "model_kind":"finite_cellular",
 "degree_range":[0,3],
 "free_ranks":[1,2,2,1],
 "boundary_maps":{
   "d1":[[enc(M1[0][0])],[enc(M1[1][0])]],
   "d2":[[enc(M2[j][i]) for i in range(2)] for j in range(2)],
   "d3":[[enc(M3[0]),enc(M3[1])]]},
 "top_closure":{"orientation":"the single 3-cell, positively oriented; its basis element is "
                "the generator of C_3 and carries the fundamental class","basis_element":0},
 "truncation_rule":None,
 "basing":{
   "module_side":"left","vector_convention":"row",
   "boundary_direction":"chains are row vectors; boundary maps act on the RIGHT, so a chain "
                        "c in C_n maps to c . d_n in C_{n-1}",
   "evaluation":"g |-> rho(g); no inverse, transpose or dual variant",
   "augmentation":"eps: C_0 -> Z sends every group element to 1; eps is the terminal map and "
                  "is NOT d_1",
   "basis_order":"degree 1 basis (e_s, e_t) matches abstract_generators; degree 2 basis "
                 "(f_1, f_2) matches the relators s^3 (st)^-2 and t^5 (st)^-2"},
 "provenance_id":{"id":"M88-CONSTR-02","source_content_sha256":hashlib.sha256(src).hexdigest()},
}
b=(json.dumps(packet,sort_keys=True,indent=2,ensure_ascii=True)+"\n").encode()
pathlib.Path("m8_8_construction_packet.json").write_bytes(b)
print(f"  {len(b)} bytes")
print(f"  packet SHA-256    : {hashlib.sha256(b).hexdigest()}")
print(f"  derivation source : {packet['provenance_id']['source_content_sha256']}")
