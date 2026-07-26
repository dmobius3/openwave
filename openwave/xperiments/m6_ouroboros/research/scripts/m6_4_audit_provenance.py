"""M6.4 ADVERSARIAL AUDIT: claim A10 (provenance: Zenodo file inventories +
printed quotes).

A10(a): live Zenodo API check (auditor's own HTTP calls, urllib, no cache):
  - record 20030162 (v1) carries exactly one file: chaoiton_theorem.lean.txt
  - record 20866581 (v2) carries exactly one file: Zenodo_20030162_Version2.docx
  - chaoiton_gf_verification.py exists on NEITHER version, although the v2
    docx's Supplementary File section claims it is attached.
A10(b): quote checks:
  - z21268405 prints "62 distinct families pass the Gelfand-Fomin stability
    test" (TEST 3, Hard Quantization) -- extracted fresh from the local docx
    with pandoc (the extract was not previously saved);
  - z20866581 section 3.2 prints "Stability is not fine-tuned: broad regions
    of parameter space pass the test".
Also: freshness check of the two SAVED fulltext extracts against fresh pandoc
runs on the archived docx files (guards every quote-based claim upstream).

Output: research/data/m6_4_audit_provenance.json
"""

import json
import subprocess
import unicodedata
import urllib.request
from pathlib import Path

HERE = Path(__file__).resolve().parent
DATA = HERE.parent / "data"
THEORY = HERE.parent.parent / "theory"

OUT = {"task": "M6.4 audit A10: provenance"}


def api(url):
    req = urllib.request.Request(url, headers={"User-Agent": "m6.4-audit"})
    with urllib.request.urlopen(req, timeout=30) as f:
        return json.loads(f.read().decode())


# ---------------------------------------------------------------- A10a: API
inv = {}
for rid in ("20030162", "20866581"):
    rec = api(f"https://zenodo.org/api/records/{rid}")
    inv[rid] = {
        "title": rec["metadata"].get("title"),
        "publication_date": rec["metadata"].get("publication_date"),
        "files": sorted(f["key"] for f in rec.get("files", [])),
    }
# version chain via the concept record search
try:
    vers = api("https://zenodo.org/api/records/20030162/versions")
    chain = [{"id": h["id"],
              "files": sorted(f["key"] for f in h.get("files", []))}
             for h in vers.get("hits", {}).get("hits", [])]
except Exception as e:  # noqa: BLE001
    chain = f"versions endpoint failed: {e}"
OUT["A10a_file_inventories"] = inv
OUT["A10a_version_chain"] = chain
allfiles = []
for v in inv.values():
    allfiles += v["files"]
if isinstance(chain, list):
    for h in chain:
        allfiles += h["files"]
OUT["A10a_gf_verification_py_anywhere"] = any(
    "chaoiton_gf_verification" in f for f in allfiles)
OUT["A10a_v1_exactly_lean"] = inv["20030162"]["files"] == [
    "chaoiton_theorem.lean.txt"]
OUT["A10a_v2_exactly_docx"] = inv["20866581"]["files"] == [
    "Zenodo_20030162_Version2.docx"]

# the v2 docx's own Supplementary File claim (from the saved extract):
v2 = (DATA / "m6_4_record_v2_fulltext.txt").read_text()
OUT["A10a_v2_supplementary_claim_present"] = (
    "available as chaoiton_gf_verification.py attached to this Zenodo record"
    in v2.replace("\n", " "))

# ------------------------------------------------------- A10b: fresh pandoc
def pandoc_plain(docx):
    r = subprocess.run(["pandoc", str(docx), "-t", "plain"],
                       capture_output=True, text=True, timeout=120)
    r.check_returncode()
    return r.stdout


def flathyph(s):
    """Flatten newlines and map unicode hyphen variants (incl. the docx's
    non-breaking hyphen U+2011) to ASCII '-'."""
    s = unicodedata.normalize("NFKC", s.replace("\n", " "))
    for h in ("‐", "‑", "‒", "–", "—"):
        s = s.replace(h, "-")
    return s


z21268405 = THEORY / ("2026-07-08_z21268405_Quantum_Measurement_Three_New_"
                      "Nuclear_Tests_v2_2_.docx")
txt3 = pandoc_plain(z21268405)
flat3 = flathyph(txt3)
needle_62 = "62 distinct families pass the Gelfand-Fomin stability test"
OUT["A10b_z21268405_62_distinct"] = needle_62 in flat3
idx = flat3.find(needle_62)
OUT["A10b_z21268405_context"] = flat3[max(0, idx - 220): idx + 120] \
    if idx >= 0 else None
OUT["A10b_z21268405_hard_quant_context"] = (
    "TEST 3 - Hard Quantization" in flat3)
OUT["A10b_z21268405_table_row"] = (
    "Hard Quantization" in flat3 and "62 stable" in flat3)

flat2 = unicodedata.normalize("NFKC", v2.replace("\n", " "))
needle_ft = ("Stability is not fine-tuned: broad regions of parameter "
             "space pass the test")
OUT["A10b_v2_not_finetuned"] = needle_ft in flat2
idx2 = flat2.find(needle_ft)
OUT["A10b_v2_not_finetuned_context"] = flat2[max(0, idx2 - 80): idx2 + 160] \
    if idx2 >= 0 else None

# ----------------------------------------- saved-extract freshness (guards)
def normtext(s):
    return " ".join(unicodedata.normalize("NFKC", s).split())


fresh_v2 = pandoc_plain(THEORY / ("2026-06-25_z20866581_Zenodo_20030162_"
                                  "Version2.docx"))
fresh_may = pandoc_plain(THEORY / ("2026-05-05_z20044392_Werbos_Chaoitons_"
                                   "Ouroboros_2025_with_far_field.docx"))
saved_may = (DATA / "m6_4_record_may_fulltext.txt").read_text()
OUT["extract_v2_matches_fresh_pandoc"] = normtext(fresh_v2) == normtext(v2)
OUT["extract_may_matches_fresh_pandoc"] = (
    normtext(fresh_may) == normtext(saved_may))

(DATA / "m6_4_audit_provenance.json").write_text(
    json.dumps(OUT, indent=2, default=str))
print(json.dumps(OUT, indent=2, default=str))
