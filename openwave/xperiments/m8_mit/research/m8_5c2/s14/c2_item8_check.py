"""S14-C2 item 8: the fallback cross-check, re-read on the C2 text. The two degraded
contracts must read the SAME at every surface that carries them, with the assignment the
C1 audit recorded: TRIGGER S items (a) and (b) at section 2 AND section 5; item (c) at
section 2 AND gate 5's row (section 5 omits (c) BY DESIGN); TRIGGER L at section 2 AND
the section 14 status classes. Armed: --selftest mutates each clause in a copy and every
mutation must go RED."""
import sys, re

def read(path):
    return open(path, encoding="utf-8").read()

def check(t):
    i=t.index("TRIGGER S"); trig=t[i:i+1800]
    s5=t[t.index("## 5. Continuation"):t.index("## 6.")]
    g5=[l for l in t.split("\n") if l.startswith("| 5 | structural identities")][0]
    s14=t[t.index("## 14. Pre-freeze"):t.index("## 15. Pins")]
    rows=[]
    def probe(name, frag, hay):
        ok = frag.lower() in hay.lower()
        rows.append((ok, name))
        return ok
    ok=True
    ok &= probe("(a) multi-seed no-symmetry @ § 2", "deterministic multi-seed search with no symmetry", trig)
    ok &= probe("(a) multi-seed fallback @ § 5", "deterministic multi-seed fallback (degraded mode)", s5)
    ok &= probe("(b) measured-nullspace @ § 2", "MEASURED-NULLSPACE", trig)
    ok &= probe("(b) rank form @ § 5", "rank_R Z", s5)
    ok &= probe("(c) energy and charge alone @ § 2", "conserved-set check runs on energy and charge alone", trig)
    ok &= probe("(c) energy and charge alone @ gate-5 row", "energy and charge alone", g5)
    ok &= probe("TRIGGER L enumeration degrade @ § 2", "enumeration degrades per (a)", trig)
    ok &= probe("TRIGGER L keeps predicted-vs-measured @ § 2", "KEEPING the predicted-versus-measured form", trig)
    ok &= probe("item 1 fallback class @ § 14", "FALLBACK-COVERED (§ 2 TRIGGER S", s14)
    ok &= probe("item 4 fallback class @ § 14", "FALLBACK-COVERED (§ 2 TRIGGER L", s14)
    for okk,name in rows: print(f"  {'OK ' if okk else 'MISS'} {name}")
    print(f"  VERDICT: {'GREEN, the same two degraded contracts at their assigned surfaces' if ok else 'RED'}")
    return ok

if "--selftest" in sys.argv:
    t=read("M8_5C2_PROTOCOL_DRAFT.md")
    muts=[("energy and charge alone","energy alone",2),
          ("MEASURED-NULLSPACE","APPROXIMATE-NULLSPACE",1),
          ("deterministic multi-seed search with no symmetry","heuristic multi-seed search with no symmetry",1)]
    bad=0
    import io, contextlib
    for o,n,cnt in muts:
        assert t.count(o)>=cnt, (o, t.count(o))
        m=t.replace(o,n)
        buf=io.StringIO()
        with contextlib.redirect_stdout(buf): fired = not check(m)
        print(f"  {'RED as required' if fired else 'FALSE-GREEN, checker broken'}: mutate {o[:40]!r}")
        bad += (not fired)
    print(f"SELFTEST: {'GREEN, every mutation fires' if bad==0 else 'RED, checker cannot fail'}")
    sys.exit(0 if bad==0 else 1)

ok=check(read("M8_5C2_PROTOCOL_DRAFT.md"))
sys.exit(0 if ok else 1)
