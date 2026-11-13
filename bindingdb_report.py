#!/usr/bin/env python3
"""
bindingdb_report.py - read-only profiling of a BindingDB TSV file -> one self-contained HTML report.

* Nothing is modified, dropped or cleaned. The script only COUNTS things and writes an HTML file.
* The report contains aggregate statistics only (counts, top-N names/IDs). No data rows are copied.
* Designed for a 16 GB RAM machine: the file is streamed in chunks and "distinct" counts use compact
  64-bit hashes stored in numpy arrays instead of Python sets.

NOTE: this code was generated with the help of an LLM (see README). Spot-check the numbers.

Usage:
    python bindingdb_report.py path/to/BindingDB_All.tsv
    python bindingdb_report.py path/to/BindingDB_All.tsv --out report.html --chunk 50000
    python bindingdb_report.py path/to/BindingDB_All.tsv --max-rows 100000     # quick trial run
"""
import argparse                      # reads the options typed after the script name
import csv                           # only used for the QUOTE_NONE constant
import html                          # html.escape() makes text safe to place inside HTML
import sys                           # sys.exit() stops the script with a message
import time                          # timing / progress output
from collections import Counter      # dict that counts things
from datetime import datetime        # timestamp for the report
from pathlib import Path             # file paths that work on Windows and Linux

import numpy as np                   # fast numeric arrays
import pandas as pd                  # chunked TSV reading

# ===================== COLUMN NAMES (as used in current BindingDB TSV files) =====================
ORG = "Target Source Organism According to Curator or DataSource"   # organism column
TNAME = "Target Name"                                                # protein name column
SMI = "Ligand SMILES"                                                # ligand structure as text
IKEY = "Ligand InChI Key"                                            # ligand structure hash-like key
SEQ = "BindingDB Target Chain Sequence 1"                            # amino-acid sequence of chain 1
NCH = "Number of Protein Chains in Target (>1 implies a multichain complex)"
AFF = ["Ki (nM)", "IC50 (nM)", "Kd (nM)", "EC50 (nM)"]               # the four affinity columns
HIST_AFF = ["Kd (nM)", "Ki (nM)", "IC50 (nM)"]                       # affinities that get a histogram
SOURCE_COLS = ["Article DOI", "PMID", "Patent Number"]               # optional: used to tell sources apart
DEFAULT_HUMAN = ["Homo sapiens", "Human", "Homo sapiens (Human)"]    # organism spellings counted as human
# =================================================================================================

PX_EDGES = np.arange(2, 13.01, 0.5)                                  # histogram bins for pKd / pKi / pIC50
SMI_EDGES = [0, 20, 40, 60, 80, 100, 150, 200, 300, 500, 1000, 10**9]    # SMILES length bins
SEQ_EDGES = [0, 100, 200, 300, 500, 750, 1000, 1500, 2000, 5000, 10**9]  # sequence length bins
SEP = "\x1f"                                                         # invisible separator used to join fields into one key


def edge_labels(edges):
    """Turn [0, 20, 40, 10**9] into ['0-20', '20-40', '40+'] for chart labels."""
    out = []
    for a, b in zip(edges[:-1], edges[1:]):                          # walk through consecutive pairs of edges
        out.append(f"{a}+" if b >= 10**9 else f"{a}-{b}")            # last bin is open-ended
    return out


def count_data_lines(path):
    """Count lines in the file (minus the header) by streaming bytes - used to detect rows pandas skipped."""
    n, last = 0, b""
    with open(path, "rb") as f:                                      # binary mode: no decoding, very fast
        while True:
            block = f.read(1 << 24)                                  # read 16 MB at a time
            if not block:                                            # empty block = end of file
                break
            n += block.count(b"\n")                                  # count newline characters
            last = block[-1:]                                        # remember the very last byte
    if last not in (b"\n", b""):                                     # file does not end with a newline
        n += 1                                                       # so the final line still counts
    return max(n - 1, 0)                                             # subtract the header line


def n_unique(arrays):
    """Number of distinct values across a list of numpy arrays."""
    return int(np.unique(np.concatenate(arrays)).size) if arrays else 0


def n_total(arrays):
    """Total number of values across a list of numpy arrays."""
    return int(sum(len(a) for a in arrays))


# ---------------------------------------------------------------- PROFILING
def profile(path, chunk, max_rows, human_labels):
    t0 = time.time()                                                 # start the clock
    header = pd.read_csv(path, sep="\t", nrows=0, quoting=csv.QUOTE_NONE).columns.tolist()  # read column names only
    need = [ORG, TNAME, SMI, IKEY, SEQ, NCH] + AFF                   # columns the report cannot work without
    miss = [c for c in need if c not in header]                      # which of them are absent
    if miss:
        sys.exit(f"Columns not found in {path}: {miss}")             # stop with a readable message
    uni = [c for c in header if "uniprot" in c.lower() and c.lower().endswith("chain 1")
           and "primary id" in c.lower()]                            # UniProt primary-ID columns for chain 1
    if not uni:
        sys.exit("No UniProt primary-ID columns (chain 1) found.")
    src_cols = [c for c in SOURCE_COLS if c in header]               # which optional source columns exist
    print(f"UniProt columns used: {uni}")
    print(f"Source columns used for the strict duplicate check: {src_cols or 'none found'}")

    file_lines = None if max_rows else count_data_lines(path)        # line count only makes sense for a full run
    if file_lines is not None:
        print(f"Data lines in file (excluding header): {file_lines:,}")

    rows = 0                                                         # rows parsed so far
    filled = np.zeros(len(header), dtype=np.int64)                   # per-column count of non-empty cells
    org_cnt = Counter()                                              # organism label -> rows
    aff = {a: Counter() for a in AFF}                                # per affinity column: blank/exact/censored/...
    hist = {a: np.zeros(len(PX_EDGES) - 1, dtype=np.int64) for a in HIST_AFF}  # histogram counts
    hist_out = Counter()                                             # values outside the histogram range
    no_aff = no_exact = human_rows = 0
    smi_blank = ikey_blank = seq_blank = uid_blank = tname_blank = 0
    h_smi, h_ikey, h_seq = [], [], []                                # lists of hash arrays (for distinct counts)
    h_a, h_b, h_c = [], [], []                                       # hash arrays for the three duplicate definitions
    smi_len = np.zeros(len(SMI_EDGES) - 1, dtype=np.int64)
    seq_len = np.zeros(len(SEQ_EDGES) - 1, dtype=np.int64)
    seq_long = seq_odd = 0
    chains, uid_rows, tname_rows = Counter(), Counter(), Counter()
    uid_seqs, uid_names = {}, {}                                     # UniProt ID -> set of sequence hashes / names

    reader = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False, quoting=csv.QUOTE_NONE,
                         encoding_errors="replace", chunksize=chunk, nrows=max_rows,
                         on_bad_lines="warn")                        # all text, no auto-NaN, no quote handling
    for ch in reader:                                                # one DataFrame of up to `chunk` rows at a time
        n = len(ch)
        rows += n
        filled += (ch != "").to_numpy().sum(axis=0)                  # non-empty cells per column

        org = ch[ORG].str.strip()                                    # organism, whitespace removed
        org_cnt.update(org.replace("", "<BLANK>").value_counts().to_dict())
        human_rows += int(org.isin(human_labels).sum())

        # ---- affinity
        any_nb = np.zeros(n, dtype=bool)                             # row has any non-blank affinity
        any_ex = np.zeros(n, dtype=bool)                             # row has any exact positive affinity
        for a in AFF:
            t = ch[a].str.strip()
            blank = (t == "").to_numpy()
            cens = t.str[:1].isin([">", "<"]).to_numpy()             # censored values start with > or <
            num = pd.to_numeric(t.str.replace(",", "", regex=False), errors="coerce").to_numpy(dtype=float)
            ok = ~np.isnan(num) & ~blank                             # parsed as a number
            pos = ok & (num > 0)                                     # ...and strictly positive
            other = ~blank & ~cens & ~ok                             # text that is neither blank, censored nor numeric
            aff[a].update(blank=int(blank.sum()), exact=int(pos.sum()), censored=int(cens.sum()),
                          other=int(other.sum()), nonpositive=int((ok & ~pos).sum()))
            any_nb |= ~blank
            any_ex |= pos
            if a in hist:
                px = 9 - np.log10(num[pos])                          # p-value = 9 - log10(value in nM)
                h = np.histogram(px, bins=PX_EDGES)[0]
                hist[a] += h
                hist_out[a] += int(len(px) - h.sum())
        no_aff += int((~any_nb).sum())
        no_exact += int((~any_ex).sum())

        # ---- ligand
        smi = ch[SMI].str.strip()
        sb = (smi == "").to_numpy()
        smi_blank += int(sb.sum())
        h_smi.append(pd.util.hash_array(smi[~sb].to_numpy(dtype=object)))   # 64-bit hashes, 8 bytes each
        smi_len += np.histogram(smi.str.len().to_numpy()[~sb], bins=SMI_EDGES)[0]
        ik = ch[IKEY].str.strip()
        ib = (ik == "").to_numpy()
        ikey_blank += int(ib.sum())
        h_ikey.append(pd.util.hash_array(ik[~ib].to_numpy(dtype=object)))

        # ---- protein
        seq = ch[SEQ].str.strip().str.upper()
        qb = (seq == "").to_numpy()
        seq_blank += int(qb.sum())
        sh = pd.util.hash_array(seq.to_numpy(dtype=object))
        h_seq.append(sh[~qb])
        sl = seq.str.len().to_numpy()
        seq_len += np.histogram(sl[~qb], bins=SEQ_EDGES)[0]
        seq_long += int((sl[~qb] > 1000).sum())
        seq_odd += int(seq[~qb].str.contains(r"[^ACDEFGHIKLMNPQRSTVWY]", regex=True).sum())

        nc = pd.to_numeric(ch[NCH].str.strip(), errors="coerce")
        lab = np.where(nc.isna(), "blank", np.where(nc >= 4, "4+", nc.fillna(0).astype(int).astype(str)))
        chains.update(pd.Series(lab).value_counts().to_dict())

        # ---- target
        uid = ch[uni[0]].str.strip()
        for c in uni[1:]:                                            # fall back to the next UniProt column if blank
            uid = uid.mask(uid == "", ch[c].str.strip())
        ub = (uid == "").to_numpy()
        uid_blank += int(ub.sum())
        uid_rows.update(uid[~ub].value_counts().to_dict())
        tn = ch[TNAME].str.strip()
        tb = (tn == "").to_numpy()
        tname_blank += int(tb.sum())
        tname_rows.update(tn[~tb].value_counts().to_dict())

        ok1 = ~ub & ~qb                                              # rows with both a UniProt ID and a sequence
        pr = pd.DataFrame({"u": uid[ok1].to_numpy(dtype=object), "h": sh[ok1]}).drop_duplicates()
        for u, h in zip(pr["u"], pr["h"]):
            uid_seqs.setdefault(u, set()).add(int(h))
        ok2 = ~ub & ~tb
        pn = pd.DataFrame({"u": uid[ok2].to_numpy(dtype=object), "n": tn[ok2].to_numpy(dtype=object)}).drop_duplicates()
        for u, nm in zip(pn["u"], pn["n"]):
            uid_names.setdefault(u, set()).add(nm)

        # ---- duplicates (only rows that have BOTH a SMILES and a UniProt ID, so blanks do not collapse together)
        keep = ~sb & ~ub
        ka = smi[keep] + SEP + uid[keep] + SEP + org[keep]           # A: ligand + target + organism
        affk = ch[AFF[0]]
        for a in AFF[1:]:
            affk = affk + SEP + ch[a]                                # join the four affinity strings
        kb = ka + SEP + affk[keep]                                   # B: A + the four affinity values
        if src_cols:
            srck = ch[src_cols[0]]
            for c in src_cols[1:]:
                srck = srck + SEP + ch[c]
            kc = kb + SEP + srck[keep]                               # C: B + source (DOI / PMID / patent)
        else:
            kc = kb
        h_a.append(pd.util.hash_array(ka.to_numpy(dtype=object)))
        h_b.append(pd.util.hash_array(kb.to_numpy(dtype=object)))
        h_c.append(pd.util.hash_array(kc.to_numpy(dtype=object)))

        print(f"  {rows:,} rows | {time.time() - t0:.0f}s", flush=True)

    dup_base = n_total(h_a)                                          # rows considered for duplicate checks
    return dict(path=str(path), header=header, rows=rows, file_lines=file_lines, filled=filled, org_cnt=org_cnt,
                aff=aff, hist=hist, hist_out=hist_out, no_aff=no_aff, no_exact=no_exact, human_rows=human_rows,
                smi_blank=smi_blank, ikey_blank=ikey_blank, seq_blank=seq_blank, uid_blank=uid_blank,
                tname_blank=tname_blank, n_smi=n_unique(h_smi), n_ikey=n_unique(h_ikey), n_seq=n_unique(h_seq),
                smi_len=smi_len, seq_len=seq_len, seq_long=seq_long, seq_odd=seq_odd, chains=chains,
                uid_rows=uid_rows, tname_rows=tname_rows, uid_seqs=uid_seqs, uid_names=uid_names,
                dup_base=dup_base, dup_a=dup_base - n_unique(h_a), dup_b=dup_base - n_unique(h_b),
                dup_c=dup_base - n_unique(h_c), src_cols=src_cols, human_labels=human_labels)


# ---------------------------------------------------------------- HTML HELPERS
def esc(x):
    return html.escape(str(x))                                       # make any text safe for HTML


def fmt(n):
    return f"{int(n):,}"                                             # 1234567 -> 1,234,567


def pct(a, b):
    return f"{100 * a / b:.1f}%" if b else "0%"


def short(n):
    n = int(n)
    return f"{n / 1e6:.1f}M" if n >= 1e6 else f"{n / 1e3:.0f}k" if n >= 1e4 else str(n)


def card(label, value, sub=""):
    return (f'<div class="card"><div class="cl">{esc(label)}</div><div class="cv">{value}</div>'
            f'<div class="cs">{esc(sub)}</div></div>')


def cards(items):
    return '<div class="cards">' + "".join(items) + "</div>"


def hbars(items, total=None, highlight=()):
    mx = max([v for _, v in items] or [1]) or 1                      # longest bar = 100% width
    out = []
    for lab, v in items:
        col = "var(--green)" if lab in highlight else "var(--blue)"
        share = f" ({100 * v / total:.1f}%)" if total else ""
        out.append(f'<div class="hb"><div class="hl" title="{esc(lab)}">{esc(lab)}</div>'
                   f'<div class="ht"><div class="hf" style="width:{100 * v / mx:.1f}%;background:{col}"></div></div>'
                   f'<div class="hv">{fmt(v)}{share}</div></div>')
    return '<div class="hbars">' + "".join(out) + "</div>"


def vbars(counts, labels, color="var(--blue)"):
    mx = max(counts) or 1
    step = 1 if len(counts) <= 12 else 2                             # label every 2nd bar when there are many
    bars, labs = [], []
    for i, (c, l) in enumerate(zip(counts, labels)):
        bars.append(f'<div class="vc" title="{esc(l)}: {fmt(c)}"><div class="vn">{short(c) if c else ""}</div>'
                    f'<div class="vb" style="height:{100 * c / mx:.1f}%;background:{color}"></div></div>')
        labs.append(f'<div class="vl">{esc(l) if i % step == 0 else ""}</div>')
    return f'<div class="vbars">{"".join(bars)}</div><div class="vlabels">{"".join(labs)}</div>'


def table(headers, rows):
    th = "".join(f"<th>{esc(h)}</th>" for h in headers)
    trs = "".join("<tr>" + "".join(f"<td>{c}</td>" for c in r) + "</tr>" for r in rows)  # cells are pre-escaped
    return f'<div class="tw"><table><thead><tr>{th}</tr></thead><tbody>{trs}</tbody></table></div>'


def section(sid, title, body, note=""):
    n = f'<p class="note">{note}</p>' if note else ""
    return f'<section id="{sid}"><h2>{title}</h2>{n}{body}</section>'


CSS = """
:root{--bg:#0f1420;--panel:#1a1f2b;--line:#2a3140;--text:#e6e9ef;--muted:#9aa4b5;--blue:#4c8dff;--green:#3fb97c}
*{box-sizing:border-box}body{margin:0;background:var(--bg);color:var(--text);font:14px/1.5 -apple-system,Segoe UI,Roboto,sans-serif}
header{padding:24px 18px;border-bottom:1px solid var(--line)}header h1{margin:0 0 6px;font-size:28px}header p{margin:0;color:var(--muted)}
nav{position:sticky;top:0;background:var(--bg);border-bottom:1px solid var(--line);padding:12px 18px;z-index:5}
nav a{color:var(--blue);text-decoration:none;margin-right:18px}
main{max-width:1440px;margin:0 auto;padding:22px 24px 60px}
section{background:var(--panel);border:1px solid var(--line);border-radius:14px;padding:24px 32px;margin-bottom:22px}
h2{margin:0 0 10px;font-size:22px}h3{margin:22px 0 8px;font-size:16px}
.note{color:var(--muted);margin:0 0 14px}
.cards{display:grid;grid-template-columns:repeat(auto-fill,minmax(210px,1fr));gap:12px;margin:14px 0}
.card{background:var(--bg);border:1px solid var(--line);border-radius:10px;padding:14px}
.cl{color:var(--muted);font-size:13px}.cv{font-size:28px;font-weight:700;margin:4px 0}.cs{color:var(--muted);font-size:12px}
.hbars{margin:10px 0}.hb{display:flex;align-items:center;gap:10px;margin:5px 0}
.hl{width:300px;flex:none;white-space:nowrap;overflow:hidden;text-overflow:ellipsis}
.ht{flex:1;background:var(--bg);border-radius:4px;height:16px}.hf{height:100%;border-radius:4px}
.hv{width:170px;flex:none;color:var(--muted);text-align:right;font-variant-numeric:tabular-nums}
.vbars{display:flex;align-items:flex-end;gap:3px;height:190px;border-bottom:1px solid var(--line);margin-top:10px}
.vc{flex:1;height:100%;display:flex;flex-direction:column;justify-content:flex-end;align-items:center}
.vb{width:100%;border-radius:3px 3px 0 0;min-height:1px}.vn{font-size:9px;color:var(--muted)}
.vlabels{display:flex;gap:3px}.vl{flex:1;text-align:center;font-size:10px;color:var(--muted);white-space:nowrap}
.tw{overflow-x:auto}table{border-collapse:collapse;width:100%;font-size:13px}
th,td{text-align:left;padding:7px 10px;border-bottom:1px solid var(--line);vertical-align:top}th{color:var(--muted);font-weight:600}
details{margin-top:12px;color:var(--muted)}summary{cursor:pointer;color:var(--blue)}
footer{max-width:1440px;margin:0 auto;padding:0 24px 40px;color:var(--muted);font-size:12px}
"""

FOOTER = (
    "<p><b>Data:</b> BindingDB (bindingdb.org). Please cite: Liu T, Hwang L, Burley SK, Nitsche CI, Southan C, "
    "Walters WP, Gilson MK. BindingDB in 2024: a FAIR knowledgebase of protein-small molecule binding data. "
    "Nucleic Acids Research 53(D1):D1633-D1644 (2025). doi:10.1093/nar/gkae1075. "
    "Data curated by BindingDB is shared under a Creative Commons Attribution license; data imported from ChEMBL "
    "is CC BY-SA 3.0 (see bindingdb.org/rwd/bind/info.jsp).</p>"
    "<p><b>About this report:</b> aggregate statistics computed on the file as provided; no rows were removed or "
    "changed and no data rows are reproduced here. The generating code was written with the help of an LLM; "
    "verify important numbers independently.</p>"
)


# ---------------------------------------------------------------- REPORT BUILDER
def build(P, title, subtitle):
    rows, header, filled = P["rows"], P["header"], P["filled"]
    ncols = len(header)
    fill_by = dict(zip(header, filled.tolist()))                     # column name -> non-empty cells
    empty = [c for c, f in zip(header, filled) if f == 0]            # columns with no data at all
    fr = filled / max(rows, 1)                                       # fill fraction per column
    low = int(((filled > 0) & (fr < 0.5)).sum())                     # non-empty columns below 50% filled
    uids = P["uid_rows"]
    multi_seq = sum(1 for s in P["uid_seqs"].values() if len(s) > 1)
    multi_name = sum(1 for s in P["uid_names"].values() if len(s) > 1)
    ch = P["chains"]
    multichain = ch.get("2", 0) + ch.get("3", 0) + ch.get("4+", 0)
    nonpos = sum(P["aff"][a]["nonpositive"] for a in AFF)
    org = P["org_cnt"]
    hl = P["human_labels"]
    human_like = [(l, n) for l, n in org.items() if "homo" in l.lower() or "human" in l.lower()]
    n_human_spellings = sum(1 for l, _ in human_like if l in hl)     # how many of your human labels actually occur
    skipped = (P["file_lines"] - rows) if P["file_lines"] is not None else None

    # 1 overview
    ov = [
        card("Rows profiled", fmt(rows), f"{ncols} columns"),
        card("Columns", fmt(ncols), f"{len(empty)} completely empty"),
        card("Columns < 50% filled", fmt(low), "excluding empty columns"),
        card("Rows matching human labels", fmt(P["human_rows"]), pct(P["human_rows"], rows) + " of rows"),
        card("Distinct ligand SMILES", fmt(P["n_smi"]), pct(P["n_smi"], rows - P["smi_blank"]) + " of non-blank"),
        card("Distinct target names", fmt(len(P["tname_rows"])), "names, not unique proteins"),
        card("Distinct UniProt IDs", fmt(len(uids)), "chain 1"),
        card("Rows with no affinity value", fmt(P["no_aff"]), pct(P["no_aff"], rows) + " (Ki, IC50, Kd, EC50 all blank)"),
        card("Rows with NO exact affinity number", fmt(P["no_exact"]), "only censored / approximate / unparseable"),
        card("Extra rows: same ligand+target+organism", fmt(P["dup_a"]), "beyond the first occurrence"),
    ]
    if skipped is not None:
        ov.append(card("Lines in file vs rows parsed", fmt(skipped),
                       f"lines not parsed ({fmt(P['file_lines'])} lines, {fmt(rows)} rows)"))
    s1 = section("overview", "1. Overview", cards(ov))

    # 2 organisms
    items = org.most_common(15)
    body = cards([card("Rows matching human labels", fmt(P["human_rows"]), " / ".join(hl)),
                  card("Rows NOT matching", fmt(rows - P["human_rows"]), "other organisms or blank"),
                  card("Blank organism rows", fmt(org.get("<BLANK>", 0)), "not counted as human"),
                  card("Distinct organism labels", fmt(len(org)), "")])
    body += "<h3>Top 15 organism labels (matched human spellings in green)</h3>" + hbars(items, rows, hl)
    top = sorted(human_like, key=lambda x: -x[1])[:25]
    body += "<h3>Every label containing 'homo' or 'human'</h3>" + table(
        ["Label", "Rows", "Counted as human?"],
        [[esc(l), fmt(n), "yes" if l in hl else "no (virus or other label)"] for l, n in top])
    s2 = section("organisms", "2. Organisms", body,
                 "Spelling variants of the same organism must be handled before filtering by species.")

    # 3 completeness
    key = [("Ligand SMILES", rows - P["smi_blank"]), ("Ligand InChI Key", rows - P["ikey_blank"]),
           ("Target sequence (chain 1)", rows - P["seq_blank"]), ("UniProt ID (chain 1)", rows - P["uid_blank"]),
           ("Target Name", rows - P["tname_blank"]), ("Organism", fill_by[ORG])] + [(a, fill_by[a]) for a in AFF]
    b = [("empty (0%)", int((fr == 0).sum())), ("0-1%", int(((fr > 0) & (fr < .01)).sum())),
         ("1-10%", int(((fr >= .01) & (fr < .1)).sum())), ("10-50%", int(((fr >= .1) & (fr < .5)).sum())),
         ("50-90%", int(((fr >= .5) & (fr < .9)).sum())), ("90-100%", int((fr >= .9).sum()))]
    body = "<h3>Key fields: rows with a value</h3>" + hbars(key, rows)
    body += "<h3>How full are the columns? (number of columns per fill band)</h3>" + hbars(b)
    s3 = section("completeness", "3. Completeness", body,
                 "A cell counts as filled if it is not an empty string (whitespace-only cells count as filled here).")

    # 4 affinity
    arows = []
    for a in AFF:
        c = P["aff"][a]
        arows.append([esc(a)] + [f"{fmt(c[k])} ({pct(c[k], rows)})" for k in ("blank", "exact", "censored", "other", "nonpositive")])
    body = cards([card("Rows with an exact Kd (>0)", fmt(P["aff"]["Kd (nM)"]["exact"]), ""),
                  card("Rows with a censored Kd", fmt(P["aff"]["Kd (nM)"]["censored"]), "> or < values"),
                  card("Rows with an exact IC50 (>0)", fmt(P["aff"]["IC50 (nM)"]["exact"]), ""),
                  card("Rows with an exact Ki (>0)", fmt(P["aff"]["Ki (nM)"]["exact"]), "")])
    body += table(["Column", "Blank", "Exact number", "Censored (> <)", "Other text", "Zero/negative"], arows)
    for a, nm in zip(HIST_AFF, ["pKd = 9 - log10(Kd nM)", "pKi = 9 - log10(Ki nM)", "pIC50 = 9 - log10(IC50 nM)"]):
        body += f"<h3>{nm} (exact values only)</h3>" + vbars(P["hist"][a].tolist(), [f"{x:.1f}" for x in PX_EDGES[:-1]])
        body += f'<p class="note">{fmt(P["hist_out"][a])} values fall outside the 2.0 to 13.0 range shown.</p>'
    s4 = section("affinity", "4. Affinity", body)

    # 5 ligand / protein
    body = cards([card("Distinct SMILES", fmt(P["n_smi"]), f"{fmt(P['smi_blank'])} blank"),
                  card("Distinct InChI Keys", fmt(P["n_ikey"]), f"{fmt(P['ikey_blank'])} blank"),
                  card("Distinct sequences", fmt(P["n_seq"]), f"{fmt(P['seq_blank'])} blank"),
                  card("Sequences > 1000 residues", fmt(P["seq_long"]), "matters for models with length limits"),
                  card("Sequences with non-standard letters", fmt(P["seq_odd"]), "outside the 20 amino acids")])
    body += "<h3>SMILES length (characters)</h3>" + vbars(P["smi_len"].tolist(), edge_labels(SMI_EDGES))
    body += "<h3>Protein sequence length (residues)</h3>" + vbars(P["seq_len"].tolist(), edge_labels(SEQ_EDGES), "var(--green)")
    body += "<h3>Number of protein chains in target</h3>" + hbars(sorted(ch.items()), rows)
    s5 = section("ligand", "5. Ligand / Protein", body,
                 "Distinct counts use 64-bit hashes of the text (collision chance is negligible at this size).")

    # 6 targets
    body = cards([card("Distinct UniProt IDs", fmt(len(uids)), f"{fmt(P['uid_blank'])} rows have none"),
                  card("UniProt IDs with >1 sequence", fmt(multi_seq), "e.g. mutants or fragments"),
                  card("UniProt IDs with >1 target name", fmt(multi_name), "naming differences")])
    body += "<h3>Top 15 target names by rows</h3>" + hbars(P["tname_rows"].most_common(15), rows)
    trs = []
    for u, n in uids.most_common(15):
        trs.append([esc(u), fmt(n), fmt(len(P["uid_seqs"].get(u, []))), esc(" | ".join(sorted(P["uid_names"].get(u, [])))[:90])])
    body += "<h3>Top 15 UniProt IDs</h3>" + table(["UniProt ID", "Rows", "Distinct sequences", "Target name(s)"], trs)
    s6 = section("targets", "6. Targets", body, "Targets are keyed by the UniProt primary ID of chain 1.")

    # 7 duplicates
    src_txt = ", ".join(P["src_cols"]) if P["src_cols"] else "none found in this file"
    drows = [
        ["A. Same ligand + UniProt target + organism", fmt(P["dup_a"]), pct(P["dup_a"], P["dup_base"]),
         "Repeated measurements of the same pair; normal in BindingDB."],
        ["B. A + identical Ki/IC50/Kd/EC50 text", fmt(P["dup_b"]), pct(P["dup_b"], P["dup_base"]),
         "Same values; may come from different papers or patents."],
        [f"C. B + identical source ({esc(src_txt)})", fmt(P["dup_c"]), pct(P["dup_c"], P["dup_base"]),
         "Most likely true repeats of one record (assay conditions are not compared)."],
    ]
    body = table(["Definition of 'duplicate'", "Extra rows", "% of rows checked", "Meaning"], drows)
    body += (f'<p class="note">Checked on {fmt(P["dup_base"])} rows that have both a SMILES and a UniProt ID. '
             'Counts are the number of rows beyond the first occurrence of each key.</p>')
    s7 = section("duplicates", "7. Duplicates", body,
                 "Three definitions, from loose to strict, because the right one depends on your question.")

    # 8 issue log
    order = {"Completeness": 1, "Structure": 2, "Validity": 3, "Duplicates": 4, "Consistency": 5}
    L = []

    def add(dim, issue, n, base, base_label, why):
        if n > 0:                                                    # only list things that were actually observed
            L.append((dim, issue, n, base, base_label, why))

    add("Completeness", "Columns that are completely empty", len(empty), ncols, "columns", "Carry no information")
    add("Completeness", "Columns less than 50% filled", low, ncols, "columns", "Sparse; check before using as features")
    add("Completeness", "Rows with no UniProt ID", P["uid_blank"], rows, "rows", "Cannot be linked to a target")
    add("Completeness", "Rows with no SMILES", P["smi_blank"], rows, "rows", "No ligand structure")
    add("Completeness", "Rows with no chain-1 sequence", P["seq_blank"], rows, "rows", "No protein sequence")
    add("Completeness", "Rows with no value in Ki/IC50/Kd/EC50", P["no_aff"], rows, "rows", "No affinity label")
    if skipped:
        add("Completeness", "Lines in the file that pandas did not parse", skipped, P["file_lines"], "lines",
            "Malformed or blank lines were skipped; inspect the file")
    add("Structure", "Multi-chain targets (>1 chain)", multichain, rows, "rows", "Complexes; only chain 1 is profiled here")
    add("Validity", "Rows with no exact affinity number", P["no_exact"], rows, "rows", "Only censored, approximate or text values")
    add("Validity", "Affinity values that are zero or negative", nonpos, rows, "rows", "Impossible for a concentration")
    add("Validity", "Sequences with non-standard letters", P["seq_odd"], rows, "rows", "Letters outside the 20 standard amino acids")
    add("Validity", "Sequences longer than 1000 residues", P["seq_long"], rows, "rows", "Relevant if a model has a length limit")
    add("Duplicates", "Extra rows: same ligand + target + organism (A)", P["dup_a"], P["dup_base"], "rows checked", "Repeated measurements")
    add("Duplicates", "Extra rows: A + identical affinity values (B)", P["dup_b"], P["dup_base"], "rows checked", "Possible repeats")
    add("Duplicates", "Extra rows: B + identical source (C)", P["dup_c"], P["dup_base"], "rows checked", "Probable repeated records")
    add("Consistency", "Rows with blank organism", org.get("<BLANK>", 0), rows, "rows", "Not counted as human")
    add("Consistency", "Human spellings found among your labels", n_human_spellings, len(hl), "labels given",
        "Several spellings of one organism")
    add("Consistency", "UniProt IDs with more than one sequence", multi_seq, max(len(uids), 1), "UniProt IDs", "Mutants may be pooled with wild type")
    add("Consistency", "UniProt IDs with more than one target name", multi_name, max(len(uids), 1), "UniProt IDs", "Naming differences")
    add("Consistency", "Distinct SMILES minus distinct InChI Keys", max(P["n_smi"] - P["n_ikey"], 0), max(P["n_smi"], 1),
        "distinct SMILES", "Same molecule written differently (if blanks are similar)")
    L.sort(key=lambda r: (order[r[0]], -r[2] / r[3]))                # group by dimension, biggest share first
    trs = [[esc(d), esc(i), fmt(n), f"{100 * n / base:.2f}%", esc(bl), esc(w)] for d, i, n, base, bl, w in L]
    s8 = section("issues", "8. Issue log",
                 table(["Dimension", "Observation", "Count", "Share", "Out of", "Why it may matter"], trs),
                 "Observations only; nothing was changed. Shares are simple proportions, not severity ratings.")

    # 9 columns
    nonempty = sorted([(c, f) for c, f in zip(header, filled) if f > 0], key=lambda x: -x[1])
    trs = [[esc(c), fmt(f), f"{100 * f / max(rows, 1):.1f}%"] for c, f in nonempty]
    body = table(["Column", "Rows filled", "% filled"], trs)
    body += (f"<details><summary>{len(empty)} completely empty columns</summary><p>"
             + esc(" | ".join(empty)) + "</p></details>")
    s9 = section("columns", "9. Columns", body)

    nav = "".join(f'<a href="#{i}">{t}</a>' for i, t in [
        ("overview", "Overview"), ("organisms", "Organisms"), ("completeness", "Completeness"), ("affinity", "Affinity"),
        ("ligand", "Ligand/Protein"), ("targets", "Targets"), ("duplicates", "Duplicates"), ("issues", "Issue log"),
        ("columns", "Columns")])
    return (f'<!doctype html><html><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">'
            f'<title>{esc(title)}</title><style>{CSS}</style></head><body>'
            f'<header><h1>{esc(title)}</h1><p>{esc(subtitle)}</p></header><nav>{nav}</nav><main>'
            f'{s1}{s2}{s3}{s4}{s5}{s6}{s7}{s8}{s9}</main><footer>{FOOTER}</footer></body></html>')


# ---------------------------------------------------------------- RUN
def main():
    ap = argparse.ArgumentParser(description="Profile a BindingDB TSV file and write a self-contained HTML report.")
    ap.add_argument("tsv", type=Path, help="path to the BindingDB TSV file (e.g. BindingDB_All.tsv)")
    ap.add_argument("--out", type=Path, default=Path("bindingdb_report.html"), help="output HTML file")
    ap.add_argument("--chunk", type=int, default=50_000, help="rows per read; lower it if you run out of RAM")
    ap.add_argument("--max-rows", type=int, default=None, help="only profile the first N rows (trial run)")
    ap.add_argument("--title", default="BindingDB data assessment", help="report title")
    ap.add_argument("--human-labels", nargs="+", default=DEFAULT_HUMAN,
                    help="organism spellings to count as human (quote labels that contain spaces)")
    args = ap.parse_args()

    if not args.tsv.exists():
        sys.exit(f"File not found: {args.tsv}")
    scope = f"first {args.max_rows:,} rows only (trial run)" if args.max_rows else "all rows in the file"
    print(f"Profiling {args.tsv.name} ({scope}), chunk = {args.chunk:,} rows")
    P = profile(args.tsv, args.chunk, args.max_rows, args.human_labels)
    stamp = datetime.now().strftime("%Y-%m-%d %H:%M")
    subtitle = f"{args.tsv.name} | {scope} | generated {stamp} | read-only, nothing removed"
    args.out.write_text(build(P, args.title, subtitle), encoding="utf-8")
    print("Saved", args.out.resolve())


if __name__ == "__main__":
    main()
