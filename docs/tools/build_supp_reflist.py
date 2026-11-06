"""Attach a self-contained numbered reference list to the IO-approval supplement.

The 137 pivotal-trial references are scoped to the supplementary table ONLY (they
are not main-manuscript Mendeley citations). This script:

  1. crosswalks each row's pivotal_trial string -> reference key(s), using the
     verification sidecar's trial names, with per-row overrides for the handful
     of trial strings that are ambiguous or composite;
  2. numbers the referenced trials in order of first appearance down the table;
  3. writes the detail CSV back with a `ref` column (the number, or "n, m" for a
     combined-trial row; blank for the 4 recent approvals with no trial paper);
  4. writes a numbered reference list (Nature-style) formatted from the verified
     Crossref/PubMed metadata in references.json.

Reference metadata is never typed by hand; it comes from references.json.

Run: python3 docs/tools/build_supp_reflist.py
"""
import csv, json, re
from pathlib import Path

root = Path(__file__).resolve().parents[1]                      # docs/
proj = root.parent
REFS = json.load(open(root / "references.json"))
meta = list(csv.DictReader(open(root / "references_meta.csv")))
detail_f = proj / "results/tables/ICI_SuppTable_io_approvals_detailed.csv"
rows = list(csv.DictReader(open(detail_f)))

# --- crosswalk: pivotal_trial string -> [keys] -------------------------------
def norm(s):
    s = re.sub(r"\(.*?\)", "", s)                               # drop parentheticals
    return re.sub(r"[^a-z0-9]", "", s.lower())

by_norm, by_full = {}, {}
for r in meta:
    by_norm.setdefault(norm(r["trial"]), []).append(r["key"])
    by_full[r["trial"].strip().lower()] = r["key"]

def split_composite(t):
    if "/" not in t:
        return [t]
    out, prefix = [], ""
    for p in (x.strip() for x in t.split("/")):
        if re.match(r"^-?\d", p) and prefix:                   # bare "-057" -> "CheckMate-057"
            out.append(f"{prefix}-{re.sub(r'^-', '', p)}")
        else:
            out.append(p)
            m = re.match(r"([A-Za-z]+)", p)
            if m:
                prefix = m.group(1)
    return out

# rows whose trial string alone is ambiguous: resolve by (trial, cancer_type)
OVERRIDE = {
    ("GARNET", "Endometrial"): ["garnet"],
    ("GARNET", "dMMR solid tumors"): ["garnetagnostic"],
    ("IND.227/KEYNOTE-483", "Mesothelioma"): ["ind227"],       # one trial w/ alias, not a split
    ("KEYNOTE-158", "Endometrial"): ["keynote158msi"],
    ("KEYNOTE-158", "Cervical"): ["keynote158cervical"],
    ("KEYNOTE-158", "TMB-high solid tumors"): ["keynote158tmb"],
    ("KEYNOTE-158/-164", "MSI-H/dMMR solid tumors"): ["keynote158msi", "keynote164"],
}

def keys_for(row):
    t = row["pivotal_trial"].strip()
    if not t:
        return []
    ov = OVERRIDE.get((t, row["cancer_type"].strip()))
    if ov:
        return ov
    if t.lower() in by_full:
        return [by_full[t.lower()]]
    keys = []
    for part in split_composite(t):
        cand = by_norm.get(norm(part), [])
        if len(cand) != 1:
            raise SystemExit(f"UNRESOLVED trial {t!r} (part {part!r}) in row {row['cancer_type']}: {cand}")
        keys.append(cand[0])
    return keys

# --- number references by order of first appearance --------------------------
order, num = [], {}
row_keys = []
for r in rows:
    ks = keys_for(r)
    row_keys.append(ks)
    for k in ks:
        if k not in num:
            order.append(k)
            num[k] = len(order)

for r, ks in zip(rows, row_keys):
    r["ref"] = ", ".join(str(num[k]) for k in ks)

# write detail CSV back with `ref` inserted after pivotal_trial
cols = list(rows[0].keys())
cols.remove("ref")
i = cols.index("pivotal_trial") + 1
cols[i:i] = ["ref"]
with open(detail_f, "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=cols, quoting=csv.QUOTE_ALL)
    w.writeheader()
    w.writerows(rows)

# --- format the numbered reference list (Nature-style) -----------------------
def initials(given):
    parts = re.split(r"[\s\-]+", given.strip())
    return " ".join(p[0].upper() + "." for p in parts if p)

def author_str(authors):
    if not authors:
        return "[author list unavailable]"
    f, g = (authors[0] + ["", ""])[:2]
    lead = f + (", " + initials(g) if g else "")
    return lead + (" et al." if len(authors) > 1 else "")

def fmt(r):
    s = f"{author_str(r['authors'])} {r['title'].rstrip('.')}. "
    if r.get("journal"):
        s += f"*{r['journal']}* "
    vp = ""
    if r.get("volume"):
        vp += f"**{r['volume']}**"
    if r.get("page"):
        vp += (", " if vp else "") + str(r["page"])
    s += vp
    s += f" ({r['year']})." if r.get("year") else "."
    if r.get("doi"):
        s += f" doi:{r['doi']}"
    return re.sub(r"\s+", " ", s).strip()

trial_of = {r["key"]: r["trial"] for r in meta}
lines = [
    "# Supplementary references: pivotal checkpoint-inhibitor approval trials",
    "",
    "Reference list for `ICI_SuppTable_io_approvals_detailed.csv` (the `ref` column "
    "keys each approved indication to its pivotal registrational trial). Scoped to "
    "the supplementary table; not part of the main-manuscript reference list. "
    "Numbered in order of first appearance in the table. Metadata resolved by DOI "
    "via Crossref/PubMed (never typed).",
    "",
]
for k in order:
    lines.append(f"{num[k]}. {fmt(REFS[k])}  <!-- {k} / {trial_of.get(k,'')} -->")

# the recent approvals with no trial paper
no_trial = [r for r, ks in zip(rows, row_keys) if not ks]
if no_trial:
    lines += ["", "## Approvals without a trial publication (regulatory source only)", ""]
    for r in no_trial:
        lines.append(f"- {r['drug']} {r['regimen']}, {r['cancer_type']} "
                     f"(FDA {r['fda_status']}, EMA {r['ema_status']}) — regulatory approval; "
                     f"pivotal-trial publication pending.")

out_f = proj / "results/tables/ICI_SuppTable_io_approvals_references.md"
out_f.write_text("\n".join(lines) + "\n")

used = set(order)
orphan = [k for k in REFS if k not in used]
print(f"rows: {len(rows)}  referenced: {sum(1 for ks in row_keys if ks)}  blank: {len(no_trial)}")
print(f"distinct references numbered: {len(order)} / {len(REFS)} in references.json")
if orphan:
    print(f"references.json keys NOT used by any row ({len(orphan)}): {sorted(orphan)}")
print(f"wrote {detail_f.name} (+ ref column) and {out_f.name}")
