"""Fetch every reference in docs/references.csv from Crossref by DOI.

Standalone metadata fetch for the immune-checkpoint-inhibitor APPROVAL-TRIAL
reference library (the pivotal registrational trials behind the FDA/EMA approval
landscape: ICI_SuppTable_io_approvals_detailed.csv). This is a separate library
from the manuscript-prose seed in manuscript/references.csv, so it is NOT coupled
to any draft: it resolves and caches every DOI it is given, cites nothing.

Reference metadata is never typed by hand: each entry is looked up by DOI and
cached in docs/references.json. Author lists come from Crossref, or from PubMed
when a row carries a PMID (the trial-group / consortium byline that Crossref
mangles is then taken from PubMed instead). PMIDs are read from the verification
sidecar docs/references_meta.csv.

Per the standing citation rule, each resolved record is checked: the DOI must be
the article itself, not a correction or erratum (Crossref type and title are
inspected; a 'journal-article' whose title is a bare page number or begins
"Correction"/"Erratum" is flagged). Flags do not type a replacement DOI; they
are reported for the user to resolve.

Run: python3 docs/tools/fetch_references.py
"""
import csv, json, re, sys, time, urllib.request
import xml.etree.ElementTree as ET
from pathlib import Path

root = Path(__file__).resolve().parents[1]                    # docs/
refs = list(csv.DictReader(open(root / "references.csv")))


def col(row, name):
    # tolerate "DOI"/"doi", "key"/"Key" headers
    for k in row:
        if k.lower() == name:
            return (row[k] or "").strip()
    return ""


# PMID per key from the verification sidecar (NA / blank -> none)
pmid_of = {}
meta_f = root / "references_meta.csv"
if meta_f.exists():
    for r in csv.DictReader(open(meta_f)):
        p = col(r, "pmid")
        if p and p.upper() != "NA":
            pmid_of[col(r, "key")] = p


def clean_title(t):
    # Crossref keeps publisher markup: "<scp>PMS2</scp> or <scp>PMS2CL</scp>? ...".
    # A one-letter small cap is a word's initial ("<scp>G</scp>\n leason") and is rejoined;
    # a longer one is a whole symbol ("<scp>PMS2</scp> or") and keeps its following space.
    t = re.sub(r"<scp>(.)</scp>\s*", r"\1", t)
    t = re.sub(r"</?scp>", "", t)
    t = re.sub(r"<[^>]+>", "", t)
    t = " ".join(t.split())
    return re.sub(r"\s+([:;,.?])", r"\1", t)


def pubmed_authors(pmid):
    url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=pubmed&id={pmid}&retmode=xml"
    tree = ET.parse(urllib.request.urlopen(url, timeout=30))
    people, groups = [], []
    for a in tree.iter("Author"):
        if a.findtext("LastName"):
            people.append((a.findtext("LastName"), a.findtext("ForeName") or a.findtext("Initials") or ""))
        elif a.findtext("CollectiveName"):
            groups.append((a.findtext("CollectiveName"), ""))
    art = tree.find(".//Article")

    def full_range(pg):                     # PubMed abbreviates: "1477-82" -> "1477-1482"
        m = re.fullmatch(r"(\d+)-(\d+)", pg or "")
        if m and len(m.group(2)) < len(m.group(1)):
            return f"{m.group(1)}-{m.group(1)[:len(m.group(1)) - len(m.group(2))]}{m.group(2)}"
        return pg or ""
    return people + groups, {"volume": art.findtext(".//JournalIssue/Volume") or "",
                             "page": full_range(art.findtext(".//Pagination/MedlinePgn"))}


def looks_like_erratum(title, typ):
    # The warning sign (per the citation rule) is a lone page number as the title,
    # or an explicit correction/erratum heading, or a Crossref correction type.
    if re.search(r"correction|erratum|author correction|retraction", typ, re.I):
        return "crossref type '%s'" % typ
    if re.fullmatch(r"\W*\d{1,5}\W*", title or ""):
        return "title is a bare page number '%s'" % title.strip()
    if re.match(r"^\s*(correction|erratum|author correction|retraction)\b", title or "", re.I):
        return "title begins '%s'" % title.strip()[:40]
    return ""


out, errors, flags = {}, [], []
for r in refs:
    k, doi = col(r, "key"), col(r, "doi")
    if not doi:
        errors.append(f"{k}: no DOI in seed")
        continue
    try:
        url = f"https://api.crossref.org/works/{doi}?mailto=mwakishasam@gmail.com"
        m = json.load(urllib.request.urlopen(url, timeout=30))["message"]
    except Exception as e:
        errors.append(f"{k}: DOI did not resolve ({doi}) - {e}")
        continue
    date = (m.get("published-print") or m.get("published-online") or m["issued"])["date-parts"][0]
    pmid = pmid_of.get(k)
    if pmid:
        authors, pm = pubmed_authors(pmid)
    else:
        authors = [(a.get("family") or a.get("name", ""), a.get("given", "")) for a in m.get("author", [])]
    title = clean_title((m.get("title") or [""])[0])
    typ = m.get("type", "")
    out[k] = {
        "doi": m["DOI"],
        "authors": authors,
        "title": title,
        "journal": (m.get("container-title") or [""])[0] or ("bioRxiv" if m["DOI"].startswith("10.1101/") else ""),
        "volume": m.get("volume", ""), "page": m.get("page", "") or m.get("article-number", ""),
        "year": date[0], "type": typ,
    }
    if pmid:                                  # fill what Crossref lacks from PubMed
        for fld in ("volume", "page"):
            out[k][fld] = out[k][fld] or pm[fld]
    why = looks_like_erratum(title, typ)
    mark = "  <-- CHECK: " + why if why else ""
    if why:
        flags.append(f"{k}: {why} (DOI {m['DOI']})")
    au = authors[0][0][:16] if authors else "?"
    print(f"{k:<22}{au:<16}{date[0]}  {title[:56]}{mark}")
    time.sleep(0.15)

json.dump(out, open(root / "references.json", "w"), indent=1, ensure_ascii=False)
print(f"\n{len(out)}/{len(refs)} references cached -> docs/references.json")
if flags:
    print(f"\n{len(flags)} possible erratum/correction DOI(s) to verify:")
    for f in flags:
        print("  " + f)
if errors:
    print(f"\n{len(errors)} unresolved row(s):")
    for e in errors:
        print("  " + e)
    sys.exit(1)
