#!/usr/bin/env python3
"""
gobi_annotate.py - find and check all members of a protein family in a genome.

One command runs the annotation pipeline from start to finish:

    genome + query proteins
      -> BLAST database (makeblastdb)
      -> tblastn search
      -> group hits into loci (one locus = one gene copy)
      -> extract each locus with flanks (samtools faidx)
      -> miniprot gene model per locus (avoids chimeras across tandem copies)
      -> checks: query coverage, frameshifts, in-frame stops, GT-AG splice sites,
         assembly gaps (N) near the locus
      -> GFF in genome coordinates (upload to the NCBI Genome Data Viewer),
         protein FASTA, summary table and a short report
      -> optional: alignment (FAMSA) + tree (IQ-TREE, UFBoot)

Usage (simplest form):

    python gobi_annotate.py --genome GCA_009733165.1 --query Q7T1K6,Q7ZZN8

  --genome  an NCBI assembly accession (GCA_/GCF_, downloaded with `datasets`)
            or a path to a genome FASTA (.fna/.fa, optionally .gz)
  --query   UniProt accessions separated by commas (downloaded automatically)
            or a path to a protein FASTA file

Every external command is printed and written to commands.log in the output
folder, so each step can be repeated by hand. Steps whose output already exists
are skipped (use --force to redo them).

Only the Python standard library is used. External tools: makeblastdb, tblastn,
samtools, miniprot; datasets (only for accessions); famsa + iqtree3 (only with --tree).
"""

import argparse
import datetime as dt
import glob
import gzip
import os
import re
import shutil
import subprocess
import sys
import urllib.request
import zipfile

__version__ = "0.1.1"

# --------------------------------------------------------------------------- #
# small helpers
# --------------------------------------------------------------------------- #

LOG = None  # open file handle for commands.log


def say(msg=""):
    print(msg, flush=True)


def die(msg):
    say(f"\nERROR: {msg}")
    sys.exit(1)


def run(cmd, stdout_path=None):
    """Run an external command, print it and log it. Abort on failure."""
    shown = " ".join(_quote(c) for c in cmd)
    if stdout_path:
        shown += f" > {_quote(stdout_path)}"
    say(f"  $ {shown}")
    if LOG:
        LOG.write(f"[{dt.datetime.now():%Y-%m-%d %H:%M:%S}] {shown}\n")
        LOG.flush()
    try:
        if stdout_path:
            with open(stdout_path, "w") as fh:
                res = subprocess.run(cmd, stdout=fh, stderr=subprocess.PIPE, text=True)
        else:
            res = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    except FileNotFoundError:
        die(f"program not found: {cmd[0]}")
    if res.returncode != 0:
        die(f"command failed ({res.returncode}): {shown}\n{(res.stderr or '')[-2000:]}")
    return res


def _quote(s):
    s = str(s)
    return f"'{s}'" if re.search(r"[\s|;&()<>]", s) else s


def done(path):
    """True if a step's output exists and is not empty."""
    return os.path.exists(path) and os.path.getsize(path) > 0


def read_fasta(path):
    """Return a list of (header, sequence)."""
    opener = gzip.open if path.endswith(".gz") else open
    recs, head, seq = [], None, []
    with opener(path, "rt") as fh:
        for line in fh:
            line = line.rstrip("\n\r")
            if line.startswith(">"):
                if head is not None:
                    recs.append((head, "".join(seq)))
                head, seq = line[1:].strip(), []
            elif line:
                seq.append(line.strip())
    if head is not None:
        recs.append((head, "".join(seq)))
    return recs


def write_fasta(path, recs, width=60):
    with open(path, "w") as fh:
        for h, s in recs:
            fh.write(f">{h}\n")
            for i in range(0, len(s), width):
                fh.write(s[i:i + width] + "\n")


COMP = str.maketrans("ACGTNacgtnRYKMSWrykmsw", "TGCANtgcanYRMKSWyrmksw")


def revcomp(s):
    return s.translate(COMP)[::-1]


_BASES = "TCAG"
_AA = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
CODON = {a + b + c: _AA[16 * i + 4 * j + k]
         for i, a in enumerate(_BASES) for j, b in enumerate(_BASES) for k, c in enumerate(_BASES)}


def translate(nt):
    nt = nt.upper()
    return "".join(CODON.get(nt[i:i + 3], "X") for i in range(0, len(nt) - 2, 3))


# --------------------------------------------------------------------------- #
# step 0: tools
# --------------------------------------------------------------------------- #

def find_iqtree():
    for name in ("iqtree3", "iqtree2", "iqtree"):
        if shutil.which(name):
            return name
    return None


def check_tools(need_datasets, need_tree):
    needed = ["makeblastdb", "tblastn", "samtools", "miniprot"]
    if need_datasets:
        needed.append("datasets")
    if need_tree:
        needed.append("famsa")
    missing = [t for t in needed if not shutil.which(t)]
    if need_tree and not find_iqtree():
        missing.append("iqtree3")
    if missing:
        die("missing tools: " + ", ".join(missing) + "\n"
            "Activate the environment first:   mamba activate gobi_annotate\n"
            "or create it once:                mamba env create -f environment.yml\n"
            "(this creates the separate environment 'gobi_annotate'; other environments are not touched)")


def tool_versions(need_tree):
    def first_line(cmd):
        try:
            r = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
            return (r.stdout.strip().splitlines() or ["?"])[0]
        except Exception:
            return "?"
    v = {
        "tblastn": first_line(["tblastn", "-version"]),
        "miniprot": first_line(["miniprot", "--version"]),
        "samtools": first_line(["samtools", "--version"]),
    }
    if need_tree:
        v["famsa"] = first_line(["famsa"])  # prints its banner
        v["iqtree"] = first_line([find_iqtree(), "--version"])
    return v


# --------------------------------------------------------------------------- #
# step 1: genome
# --------------------------------------------------------------------------- #

ACC_RE = re.compile(r"^GC[AF]_\d{9}\.\d+$")


def prepare_genome(genome_arg, cache_root, force):
    """Return (name, path to uncompressed genome FASTA, cache dir)."""
    if os.path.exists(genome_arg):
        base = os.path.basename(genome_arg)
        name = re.sub(r"(\.gz)?$", "", base)
        name = re.sub(r"\.(fna|fa|fasta|fas)$", "", name)
        name = re.sub(r"[^A-Za-z0-9_.-]", "_", name)  # cache/BLAST names without spaces
        cdir = os.path.join(cache_root, name)
        os.makedirs(cdir, exist_ok=True)
        if genome_arg.endswith(".gz"):
            fna = os.path.join(cdir, name + ".fna")
            if force or not done(fna):
                say(f"  unpacking {genome_arg}")
                with gzip.open(genome_arg, "rb") as src, open(fna, "wb") as dst:
                    shutil.copyfileobj(src, dst)
            return name, fna, cdir
        return name, os.path.abspath(genome_arg), cdir

    if not ACC_RE.match(genome_arg):
        die(f"--genome is neither an existing file nor an NCBI accession: {genome_arg}")
    name = genome_arg
    cdir = os.path.join(cache_root, name)
    os.makedirs(cdir, exist_ok=True)
    fna = os.path.join(cdir, name + ".fna")
    if not force and done(fna):
        say(f"  genome already downloaded: {fna}")
        return name, fna, cdir
    zpath = os.path.join(cdir, name + ".zip")
    run(["datasets", "download", "genome", "accession", name,
         "--include", "genome", "--filename", zpath])
    with zipfile.ZipFile(zpath) as z:
        members = [m for m in z.namelist() if m.endswith(".fna")]
        if not members:
            die(f"no .fna file found in {zpath}")
        say(f"  extracting {members[0]}")
        with z.open(members[0]) as src, open(fna, "wb") as dst:
            shutil.copyfileobj(src, dst)
    os.remove(zpath)
    return name, fna, cdir


def read_fai(fna):
    lens = {}
    with open(fna + ".fai") as fh:
        for line in fh:
            f = line.split("\t")
            lens[f[0]] = int(f[1])
    return lens


# --------------------------------------------------------------------------- #
# step 2: queries
# --------------------------------------------------------------------------- #

def prepare_queries(query_arg, out_fa):
    if os.path.exists(query_arg):
        recs = read_fasta(query_arg)
    else:
        accs = [a.strip() for a in query_arg.split(",") if a.strip()]
        recs = []
        for acc in accs:
            url = f"https://rest.uniprot.org/uniprotkb/{acc}.fasta"
            say(f"  downloading {url}")
            if LOG:
                LOG.write(f"[{dt.datetime.now():%Y-%m-%d %H:%M:%S}] GET {url}\n")
            try:
                req = urllib.request.Request(url, headers={"User-Agent": "gobi_annotate"})
                with urllib.request.urlopen(req, timeout=60) as r:
                    txt = r.read().decode()
            except Exception as e:
                die(f"could not download {acc} from UniProt: {e}")
            if not txt.startswith(">"):
                die(f"UniProt returned no sequence for {acc}")
            lines = txt.strip().splitlines()
            recs.append((lines[0][1:], "".join(lines[1:])))
    if not recs:
        die("no query sequences")
    recs = [(h, s.replace("*", "").upper()) for h, s in recs]
    write_fasta(out_fa, recs)
    return {h.split()[0]: len(s) for h, s in recs}


# --------------------------------------------------------------------------- #
# step 3-4: BLAST
# --------------------------------------------------------------------------- #

def build_blastdb(fna, cdir, name, force):
    db = os.path.join(cdir, "blastdb", name)
    if " " in db:
        say("  WARNING: the BLAST database path contains a space; BLAST may fail.\n"
            "           Use --cache with a path without spaces.")
    os.makedirs(os.path.dirname(db), exist_ok=True)
    if not force and (glob.glob(db + "*.nsq") or glob.glob(db + "*.nal")):
        say(f"  BLAST database exists: {db}")
        return db
    src = fna
    if " " in fna:
        # makeblastdb reads -in as a space-separated list of files, so a path with
        # spaces (e.g. iCloud "Mobile Documents") breaks it. Use a link without spaces.
        src = os.path.join(cdir, name + "_link.fna")
        if os.path.lexists(src):
            os.remove(src)
        os.symlink(fna, src)
        say(f"  genome path contains spaces -> using link {src}")
    run(["makeblastdb", "-in", src, "-dbtype", "nucl", "-parse_seqids", "-out", db])
    return db


HIT_COLS = ["qseqid", "saccver", "pident", "length", "qstart", "qend",
            "sstart", "send", "evalue", "bitscore"]


def run_tblastn(queries, db, out_tsv, evalue, threads, force):
    if not force and done(out_tsv):
        say(f"  tblastn output exists: {out_tsv}")
        return
    run(["tblastn", "-query", queries, "-db", db, "-evalue", str(evalue),
         "-num_threads", str(threads), "-outfmt", "6 " + " ".join(HIT_COLS),
         "-out", out_tsv])


def read_hits(tsv):
    hits = []
    if not os.path.exists(tsv):
        return hits
    with open(tsv) as fh:
        for line in fh:
            if not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            h = dict(zip(HIT_COLS, f))
            for k in ("qstart", "qend", "sstart", "send", "length"):
                h[k] = int(h[k])
            for k in ("pident", "evalue", "bitscore"):
                h[k] = float(h[k])
            h["strand"] = "+" if h["sstart"] <= h["send"] else "-"
            h["gstart"] = min(h["sstart"], h["send"])
            h["gend"] = max(h["sstart"], h["send"])
            hits.append(h)
    return hits


# --------------------------------------------------------------------------- #
# step 5: group hits into loci
# --------------------------------------------------------------------------- #

def cluster_loci(hits, max_intron, tol=20):
    """
    One locus = one gene copy. Hits on the same contig and strand are joined
    while (a) the gap to the previous hit is <= max_intron and (b) the query
    position keeps moving in the direction expected for that strand. If the
    query "starts again" (e.g. exon 1 appears after the last exon), a new copy
    begins - this separates tandem copies that lie close together.
    """
    loci = []
    hits = sorted(hits, key=lambda h: (h["saccver"], h["strand"], h["gstart"]))
    cur = None
    for h in hits:
        new = cur is None or h["saccver"] != cur["contig"] or h["strand"] != cur["strand"] \
            or h["gstart"] - cur["end"] > max_intron
        if not new:
            last = cur["last"].get(h["qseqid"])
            if last is not None:
                if h["strand"] == "+" and h["qstart"] < last["qstart"] - tol:
                    new = True
                if h["strand"] == "-" and h["qstart"] > last["qstart"] + tol:
                    new = True
        if new:
            cur = {"contig": h["saccver"], "strand": h["strand"], "start": h["gstart"],
                   "end": h["gend"], "hits": [], "last": {}}
            loci.append(cur)
        cur["hits"].append(h)
        cur["end"] = max(cur["end"], h["gend"])
        cur["last"][h["qseqid"]] = h
    loci.sort(key=lambda L: (L["contig"], L["start"]))
    width = max(2, len(str(len(loci))))
    for i, L in enumerate(loci, 1):
        L["id"] = f"L{i:0{width}d}"
        L["best_evalue"] = min(h["evalue"] for h in L["hits"])
        L["queries"] = sorted({h["qseqid"] for h in L["hits"]})
    return loci


def set_windows(loci, contig_len, flank):
    """Window = locus +/- flank, but never reaching into a neighbouring locus."""
    by_contig = {}
    for L in loci:
        by_contig.setdefault(L["contig"], []).append(L)
    for contig, Ls in by_contig.items():
        Ls.sort(key=lambda L: L["start"])
        for i, L in enumerate(Ls):
            lo = max(1, L["start"] - flank)
            hi = min(contig_len[contig], L["end"] + flank)
            if i > 0:
                lo = max(lo, Ls[i - 1]["end"] + 1)
            if i < len(Ls) - 1:
                hi = min(hi, Ls[i + 1]["start"] - 1)
            L["win_start"], L["win_end"] = lo, hi


# --------------------------------------------------------------------------- #
# step 6-7: miniprot per locus + checks
# --------------------------------------------------------------------------- #

def parse_attrs(s):
    d = {}
    for part in s.strip().split(";"):
        if "=" in part:
            k, v = part.split("=", 1)
            d[k] = v
    return d


def read_miniprot_gff(path):
    """Return list of models: dict(mrna fields, cds list, has_stop_feature)."""
    models, by_id = [], {}
    pending = {}  # ##PAF / ##STA lines come right before the model they belong to
    with open(path) as fh:
        for line in fh:
            if line.startswith("##PAF"):
                tags = dict(t.split(":i:") for t in line.split("\t") if ":i:" in t)
                pending = {"fs": int(tags.get("fs", 0)), "st": int(tags.get("st", 0))}
                continue
            if line.startswith("##STA"):
                pending["sta"] = line.rstrip("\n").split("\t")[1]
                continue
            if line.startswith("#") or not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9:
                continue
            a = parse_attrs(f[8])
            if f[2] == "mRNA":
                t = a.get("Target", "").split()
                m = {"seqid": f[0], "start": int(f[3]), "end": int(f[4]),
                     "score": float(f[5]) if f[5] != "." else 0.0, "strand": f[6],
                     "id": a.get("ID"), "identity": float(a.get("Identity", 0)),
                     "positive": float(a.get("Positive", 0)),
                     "frameshift": int(a.get("Frameshift", 0)),
                     "stopcodon_attr": int(a.get("StopCodon", 0)),
                     "target": t[0] if t else "?",
                     "tstart": int(t[1]) if len(t) > 2 else 0,
                     "tend": int(t[2]) if len(t) > 2 else 0,
                     "cds": [], "has_stop": False, "lines": [f],
                     "sta": pending.get("sta"), "paf_fs": pending.get("fs", 0),
                     "paf_st": pending.get("st", 0)}
                pending = {}
                models.append(m)
                by_id[m["id"]] = m
            else:
                m = by_id.get(a.get("Parent"))
                if m is None:
                    continue
                m["lines"].append(f)
                if f[2] == "CDS":
                    m["cds"].append((int(f[3]), int(f[4])))
                elif f[2] == "stop_codon":
                    m["has_stop"] = True
    return models


def check_model(m, winseq, qlen, min_intron=20):
    """Translate the model and run the by-eye checks automatically."""
    cds = sorted(m["cds"])
    nt = "".join(winseq[a - 1:b] for a, b in cds)
    introns = [(cds[i][1] + 1, cds[i + 1][0] - 1) for i in range(len(cds) - 1)]
    real_introns = [(a, b) for a, b in introns if b - a + 1 >= min_intron]
    short_gaps = len(introns) - len(real_introns)  # 1-19 nt gaps = frameshift fixes
    if m["strand"] == "-":
        nt = revcomp(nt)
    # miniprot's own translation (--trans) skips frameshifts correctly; use it when present
    prot = m["sta"] if m.get("sta") else translate(nt)
    ends_with_stop = prot.endswith("*")
    prot_core = prot.rstrip("*")
    internal_stops = prot_core.count("*")

    bad_splice = 0
    for a, b in real_introns:
        s = winseq[a - 1:b].upper()
        if m["strand"] == "-":
            s = revcomp(s)
        if not ((s.startswith("GT") or s.startswith("GC")) and s.endswith("AG")):
            bad_splice += 1

    span = winseq[m["start"] - 1:m["end"]].upper()
    cov = (m["tend"] - m["tstart"] + 1) / qlen if qlen else 0
    return {
        "protein": prot_core,
        "aa_len": len(prot_core),
        "coverage": cov,
        "n_exons": len(real_introns) + 1,
        "introns": len(real_introns),
        "bad_splice": bad_splice,
        "frameshift": max(m["frameshift"], m["paf_fs"], short_gaps),
        "internal_stops": max(internal_stops, m["stopcodon_attr"], m["paf_st"]),
        "has_stop": m["has_stop"] or ends_with_stop,
        "starts_M": prot_core.startswith("M"),
        "n_in_model": span.count("N"),
        "max_intron": max([b - a + 1 for a, b in real_introns], default=0),
        "cys": prot_core.count("C"),
    }


def verdict(c, n_in_window, min_cov):
    problems = []
    if c["frameshift"]:
        problems.append(f"frameshift x{c['frameshift']}")
    if c["internal_stops"]:
        problems.append(f"in-frame stop x{c['internal_stops']}")
    if c["bad_splice"]:
        problems.append(f"non-GT-AG intron x{c['bad_splice']}")
    if c["coverage"] < min_cov:
        problems.append(f"covers {c['coverage']:.0%} of query")
    if not c["starts_M"]:
        problems.append("no start M")
    if not c["has_stop"]:
        problems.append("no stop codon")
    if n_in_window:
        problems.append(f"{n_in_window} N (assembly gap) in window")
    if c["frameshift"] or c["internal_stops"]:
        v = "suspicious"   # possible pseudogene or sequencing/assembly error
    elif c["coverage"] >= min_cov and not c["bad_splice"] and c["starts_M"] and c["has_stop"]:
        v = "complete"
    else:
        v = "partial"
    return v, "; ".join(problems) if problems else "-"


def process_loci(loci, fna, queries_fa, qlens, out, threads, force, min_cov):
    ldir = os.path.join(out, "loci")
    os.makedirs(ldir, exist_ok=True)
    all_gff, best_gff, proteins = [], [], []
    for L in loci:
        region = f"{L['contig']}:{L['win_start']}-{L['win_end']}"
        wfa = os.path.join(ldir, f"{L['id']}.fa")
        wgff = os.path.join(ldir, f"{L['id']}_miniprot.gff")
        if force or not done(wfa):
            run(["samtools", "faidx", fna, region], stdout_path=wfa)
        if force or not done(wgff):
            run(["miniprot", f"-t{threads}", "-N", "20", "--outs=0.3", "-p", "0.3",
                 "--gff", "--trans", wfa, queries_fa], stdout_path=wgff)
        winseq = read_fasta(wfa)[0][1]
        L["n_in_window"] = winseq.upper().count("N")
        offset = L["win_start"] - 1
        models = read_miniprot_gff(wgff)
        # keep models that overlap the tblastn hits of this locus and are on its strand
        core = (L["start"] - offset, L["end"] - offset)
        models = [m for m in models if m["strand"] == L["strand"]
                  and m["start"] <= core[1] and m["end"] >= core[0]]
        for m in models:
            for f in m["lines"]:
                g = list(f)
                g[0] = L["contig"]
                g[3] = str(int(g[3]) + offset)
                g[4] = str(int(g[4]) + offset)
                g[8] = re.sub(r"(ID|Parent)=", lambda x: x.group(0) + L["id"] + "_", g[8])
                all_gff.append(g)
        if not models:
            L.update(verdict="no_model", notes="miniprot found no model", best=None)
            continue
        best = max(models, key=lambda m: m["score"])
        c = check_model(best, winseq, qlens.get(best["target"], 0))
        v, notes = verdict(c, L["n_in_window"], min_cov)
        L.update(best=best, check=c, verdict=v, notes=notes,
                 model_start=best["start"] + offset, model_end=best["end"] + offset)
        for f in best["lines"]:
            g = list(f)
            g[0] = L["contig"]
            g[3] = str(int(g[3]) + offset)
            g[4] = str(int(g[4]) + offset)
            g[8] = re.sub(r"(ID|Parent)=", lambda x: x.group(0) + L["id"] + "_", g[8])
            g[8] += f";locus={L['id']};verdict={v}"
            best_gff.append(g)
        short_q = best["target"].split("|")[-1]
        proteins.append((f"{L['id']}|{L['contig']}:{L['model_start']}-{L['model_end']}({L['strand']})"
                         f"|{v}|best_query={short_q}", c["protein"]))

    def write_gff(path, rows):
        with open(path, "w") as fh:
            fh.write("##gff-version 3\n")
            for g in rows:
                fh.write("\t".join(g) + "\n")
    write_gff(os.path.join(out, "models_all.gff"), all_gff)
    write_gff(os.path.join(out, "models_best.gff"), best_gff)
    write_fasta(os.path.join(out, "proteins.faa"), proteins)
    return proteins


# --------------------------------------------------------------------------- #
# step 8: tables and report
# --------------------------------------------------------------------------- #

SUMMARY_COLS = ["locus", "contig", "strand", "hits_start", "hits_end", "n_hits", "best_evalue",
                "model_start", "model_end", "best_query", "identity", "query_from", "query_to",
                "coverage", "aa_len", "exons", "max_intron", "cys", "frameshift",
                "inframe_stops", "non_gtag_introns", "N_in_window", "verdict", "notes"]


def write_tables(loci, out):
    with open(os.path.join(out, "hits_by_locus.tsv"), "w") as fh:
        fh.write("locus\t" + "\t".join(HIT_COLS) + "\n")
        for L in loci:
            for h in L["hits"]:
                fh.write(L["id"] + "\t" + "\t".join(str(h[k]) for k in HIT_COLS) + "\n")
    rows = []
    for L in loci:
        b, c = L.get("best"), L.get("check", {})
        r = {"locus": L["id"], "contig": L["contig"], "strand": L["strand"],
             "hits_start": L["start"], "hits_end": L["end"], "n_hits": len(L["hits"]),
             "best_evalue": f"{L['best_evalue']:.1e}", "verdict": L["verdict"],
             "notes": L["notes"], "N_in_window": L.get("n_in_window", "")}
        if b:
            r.update(model_start=L["model_start"], model_end=L["model_end"],
                     best_query=b["target"], identity=f"{b['identity']:.3f}",
                     query_from=b["tstart"], query_to=b["tend"],
                     coverage=f"{c['coverage']:.2f}", aa_len=c["aa_len"], exons=c["n_exons"],
                     max_intron=c["max_intron"], cys=c["cys"], frameshift=c["frameshift"],
                     inframe_stops=c["internal_stops"], non_gtag_introns=c["bad_splice"])
        rows.append(r)
    with open(os.path.join(out, "summary.tsv"), "w") as fh:
        fh.write("\t".join(SUMMARY_COLS) + "\n")
        for r in rows:
            fh.write("\t".join(str(r.get(k, "")) for k in SUMMARY_COLS) + "\n")
    return rows


def write_report(out, args, gname, qlens, hits, loci, rows, versions, tree_info):
    counts = {}
    for r in rows:
        counts[r["verdict"]] = counts.get(r["verdict"], 0) + 1
    lines = [
        f"# gobi_annotate report: {gname}",
        "",
        f"- Date: {dt.datetime.now():%Y-%m-%d %H:%M}",
        f"- Pipeline version: {__version__}",
        f"- Command: `{' '.join(_quote(a) for a in sys.argv)}`",
        f"- Genome: {args.genome}",
        f"- Queries ({len(qlens)}): " + ", ".join(f"{q} ({n} aa)" for q, n in qlens.items()),
        f"- Settings: E-value {args.evalue}, max intron {args.max_intron} bp, "
        f"flank {args.flank} bp, complete if coverage >= {args.min_cov:.0%}",
        "- Tools: " + "; ".join(f"{k}: {v}" for k, v in versions.items()),
        "",
        "## Result",
        "",
        f"- tblastn hits: {len(hits)}",
        f"- loci (gene copies): {len(loci)}",
        "- verdicts: " + (", ".join(f"{k} {v}" for k, v in sorted(counts.items())) or "none"),
        "",
        "| locus | position | strand | best query | identity | query aa | aa | exons | verdict | notes |",
        "|---|---|---|---|---|---|---|---|---|---|",
    ]
    for r in rows:
        pos = f"{r['contig']}:{r.get('model_start') or r['hits_start']:,}-{r.get('model_end') or r['hits_end']:,}"
        q = str(r.get("best_query", "-")).split("|")[-1]
        qa = f"{r.get('query_from', '')}-{r.get('query_to', '')}" if r.get("best_query") else "-"
        lines.append(f"| {r['locus']} | {pos} | {r['strand']} | {q} | {r.get('identity', '-')} | "
                     f"{qa} | {r.get('aa_len', '-')} | {r.get('exons', '-')} | **{r['verdict']}** | {r['notes']} |")
    lines += [
        "",
        "Verdicts: **complete** = covers the query, starts with M, ends with a stop, all introns GT-AG, "
        "no frameshift; **partial** = part of the query missing or not clean at the ends; "
        "**suspicious** = frameshift or in-frame stop (pseudogene, or sequencing/assembly error); "
        "**no_model** = tblastn hit but miniprot could not build a model.",
        "",
        "These are automatic checks. Look at the loci by eye before calling a copy a gene "
        "(\"not found\" is not \"absent\", and one hit is not one gene).",
        "",
        "## Files",
        "",
        "- `summary.tsv`: one row per locus (all checks)",
        "- `hits_by_locus.tsv`: every tblastn hit with its locus",
        "- `models_best.gff`: best miniprot model per locus, genome coordinates",
        "- `models_all.gff`: all miniprot models (all queries), genome coordinates",
        "- `proteins.faa`: translated best model per locus",
        "- `loci/`: extracted window and raw miniprot output per locus",
        "- `commands.log`: every command that was run, with time stamps",
    ]
    if tree_info:
        lines += [f"- `tree/`: {tree_info}"]
    lines += [
        "",
        "## Look at it by eye (NCBI Genome Data Viewer)",
        "",
        "Open the assembly in the NCBI Genome Data Viewer, then Tracks > Configure Tracks > "
        "Custom Data > Local File > upload `models_best.gff`. The contig names in the GFF are the "
        "sequence names of the genome FASTA (for NCBI assemblies the GenBank accession, e.g. CM019148.1).",
    ]
    with open(os.path.join(out, "report.md"), "w") as fh:
        fh.write("\n".join(lines) + "\n")


# --------------------------------------------------------------------------- #
# step 9: optional tree
# --------------------------------------------------------------------------- #

def safe_name(s):
    return re.sub(r"[^A-Za-z0-9_.-]", "_", s)[:60]


def make_tree(out, queries_fa, proteins, refs, threads, bootstrap, model, force):
    tdir = os.path.join(out, "tree")
    os.makedirs(tdir, exist_ok=True)
    recs, seen = [], set()

    def add(h, s):
        n = safe_name(h.split()[0])
        while n in seen:
            n += "_x"
        seen.add(n)
        recs.append((n, s.replace("*", "")))
    for h, s in read_fasta(queries_fa):
        add(h, s)
    if refs:
        for h, s in read_fasta(refs):
            add(h, s)
    for h, s in proteins:
        if len(s) >= 30:
            add(h.split("|")[0] + "_" + h.split("|")[2], s)
    if len(recs) < 4:
        say("  fewer than 4 sequences - no tree (IQ-TREE needs at least 4)")
        return "skipped (fewer than 4 sequences)"
    inp = os.path.join(tdir, "tree_input.faa")
    aln = os.path.join(tdir, "alignment.faa")
    write_fasta(inp, recs)
    if force or not done(aln):
        run(["famsa", "-t", str(threads), inp, aln])
    prefix = os.path.join(tdir, "tree")
    if force or not done(prefix + ".treefile"):
        cmd = [find_iqtree(), "-s", aln, "-m", model, "-B", str(bootstrap),
               "-T", str(max(1, min(threads, os.cpu_count() or 1))), "--prefix", prefix, "-redo"]
        run(cmd)
    return (f"{len(recs)} sequences; `alignment.faa` (FAMSA), `tree.treefile` "
            f"(IQ-TREE, model {model}, {bootstrap} UFBoot) - open in iTOL")


# --------------------------------------------------------------------------- #
# main
# --------------------------------------------------------------------------- #

class _Help(argparse.ArgumentDefaultsHelpFormatter):
    """Show defaults, but not for options without one or for on/off switches."""
    def _get_help_string(self, action):
        if action.default in (None, False) or action.default is argparse.SUPPRESS:
            return action.help
        return super()._get_help_string(action)


def main():
    p = argparse.ArgumentParser(
        description="Find and check all copies of a protein family in a genome "
                    "(tblastn + miniprot + automatic checks).",
        formatter_class=_Help)
    p.add_argument("--genome", required=True,
                   help="NCBI assembly accession (e.g. GCA_009733165.1) or genome FASTA file")
    p.add_argument("--query", required=True,
                   help="UniProt accessions, comma-separated (e.g. Q7T1K6,Q7ZZN8), or a protein FASTA")
    p.add_argument("--out", help="output folder (default: runs/<genome>_<date>)")
    p.add_argument("--tree", action="store_true", help="also build alignment (FAMSA) and tree (IQ-TREE)")
    p.add_argument("--refs", help="extra reference proteins (FASTA) for the tree, e.g. outgroups")
    p.add_argument("--threads", type=int, default=4, help="CPU threads")
    p.add_argument("--evalue", type=float, default=1e-3, help="tblastn E-value cut-off")
    p.add_argument("--max-intron", type=int, default=20000,
                   help="largest gap (bp) between two hits of the same gene copy")
    p.add_argument("--flank", type=int, default=10000,
                   help="bp added on each side of a locus before running miniprot")
    p.add_argument("--min-cov", type=float, default=0.9,
                   help="query coverage needed to call a copy complete")
    p.add_argument("--model", default="MFP", help="IQ-TREE model (MFP = ModelFinder)")
    p.add_argument("--bootstrap", type=int, default=1000, help="UFBoot replicates")
    p.add_argument("--cache", default=os.path.join(os.path.expanduser("~"), "gobi_cache"),
                   help="where genomes and BLAST databases are kept (path without spaces)")
    p.add_argument("--force", action="store_true", help="redo all steps even if output exists")
    p.add_argument("--version", action="version", version=__version__)
    args = p.parse_args()

    global LOG
    is_file = os.path.exists(args.genome)
    check_tools(need_datasets=not is_file, need_tree=args.tree)

    gtag = re.sub(r"\.(fna|fa|fasta|fas)(\.gz)?$", "", os.path.basename(args.genome))
    gtag = re.sub(r"[^A-Za-z0-9_.-]", "_", gtag)
    out = args.out or os.path.join("runs", f"{gtag}_{dt.date.today():%Y-%m-%d}")
    os.makedirs(out, exist_ok=True)
    LOG = open(os.path.join(out, "commands.log"), "a")
    LOG.write(f"\n# run started {dt.datetime.now():%Y-%m-%d %H:%M:%S}: "
              f"{' '.join(_quote(a) for a in sys.argv)}\n")
    say(f"gobi_annotate {__version__}  ->  output folder: {out}\n")

    say("[1/7] genome")
    gname, fna, cdir = prepare_genome(args.genome, args.cache, args.force)
    if args.force or not done(fna + ".fai"):
        run(["samtools", "faidx", fna])
    contig_len = read_fai(fna)
    say(f"  {len(contig_len):,} sequences, {sum(contig_len.values()):,} bp")

    say("[2/7] queries")
    qfa = os.path.join(out, "queries.faa")
    qlens = prepare_queries(args.query, qfa)
    for q, n in qlens.items():
        say(f"  {q}  {n} aa")

    say("[3/7] BLAST database")
    db = build_blastdb(fna, cdir, gname, args.force)

    say("[4/7] tblastn")
    hits_tsv = os.path.join(out, "tblastn_hits.tsv")
    run_tblastn(qfa, db, hits_tsv, args.evalue, args.threads, args.force)
    hits = read_hits(hits_tsv)
    say(f"  {len(hits)} hits")
    unknown = {h["saccver"] for h in hits} - set(contig_len)
    if unknown:
        die("tblastn reports sequence names that are not in the genome index: "
            + ", ".join(sorted(unknown)[:5]))

    say("[5/7] group hits into loci")
    loci = cluster_loci(hits, args.max_intron)
    set_windows(loci, contig_len, args.flank)
    say(f"  {len(loci)} loci")

    say("[6/7] miniprot per locus + checks")
    proteins = process_loci(loci, fna, qfa, qlens, out, args.threads, args.force, args.min_cov) \
        if loci else []
    rows = write_tables(loci, out)

    tree_info = None
    if args.tree:
        say("[7/7] alignment + tree")
        tree_info = make_tree(out, qfa, proteins, args.refs, args.threads,
                              args.bootstrap, args.model, args.force)
    else:
        say("[7/7] tree skipped (add --tree)")

    write_report(out, args, gname, qlens, hits, loci, rows, tool_versions(args.tree), tree_info)

    say("\nResult:")
    for r in rows:
        q = str(r.get("best_query", "-")).split("|")[-1]
        say(f"  {r['locus']}  {r['contig']}:{r['hits_start']}-{r['hits_end']} ({r['strand']})  "
            f"{q:<14} {r['verdict']:<10} {r['notes']}")
    say(f"\nDone. Open {os.path.join(out, 'report.md')}")
    LOG.close()


if __name__ == "__main__":
    main()
