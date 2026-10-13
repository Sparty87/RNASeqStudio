#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
RNAseqStudio — complete RNA-seq analysis with a graphical interface.

Steps, each usable on its own:
  1. FASTQ  -> (fastp trimming) -> HISAT2 / STAR alignment -> sorted, indexed BAM files
  2. BAM    -> gene-level counts with featureCounts (+ TPM and gene annotation from the GTF)
  3. Count matrix -> choice of samples and controls -> differential expression (PyDESeq2)
  4. DE results -> pathway analysis (hypergeometric ORA and preranked GSEA)

Self-installation (no administrator rights, no virtual environment):
  - Python dependencies: installed for the current user in ~/.local on first launch;
  - external tools (fastp, HISAT2, STAR, subread/featureCounts, samtools): Bioconda environment
    in ~/.local/share/rnaseqstudio/tools, installed from the "Environment" page or with --install-tools.
  Steps 3 and 4 do not require any external tool.

Author: Serena Spartano — Università Cattolica del Sacro Cuore, Rome, Italy
License: MIT

Usage:
  python3 rnaseq_studio.py                  start the GUI
  python3 rnaseq_studio.py --install-tools  install the external tools from the command line (no display needed)
"""
from __future__ import annotations

import os
import sys
import subprocess
import importlib.util
import platform
import shutil
import sysconfig
import urllib.request
from pathlib import Path

APP_NAME = "RNAseqStudio"
APP_VERSION = "1.0.0"
__version__ = APP_VERSION
__author__ = "Serena Spartano"
APP_HOME = Path(os.environ.get("RNASEQSTUDIO_HOME",
                               str(Path.home() / ".local" / "share" / "rnaseqstudio"))).expanduser()
LOCAL_BIN = Path.home() / ".local" / "bin"
TOOLS_ENV = APP_HOME / "tools"
MAMBA_BIN = APP_HOME / "bin" / "micromamba"
GENESET_CACHE = APP_HOME / "genesets"

PY_DEPS = {  # importable module -> pip requirement
    "PySide6": "PySide6>=6.5",
    "pandas": "pandas>=1.5",
    "numpy": "numpy>=1.23",
    "scipy": "scipy>=1.9",
    "matplotlib": "matplotlib>=3.6",
    "pydeseq2": "pydeseq2>=0.4",
    "gseapy": "gseapy>=1.0",
}
CONDA_PKGS = ["fastp", "hisat2", "star", "subread", "samtools"]
EXT_TOOLS = ["fastp", "hisat2", "hisat2-build", "STAR", "featureCounts", "samtools"]


# =============================================================================
# 0. External tools and bootstrap (standard library only)
# =============================================================================
def tool_path(name: str):
    p = TOOLS_ENV / "bin" / name
    if p.exists():
        return str(p)
    return shutil.which(name)


def tool_env() -> dict:
    env = os.environ.copy()
    env["MAMBA_ROOT_PREFIX"] = str(APP_HOME / "mamba_root")
    if (TOOLS_ENV / "bin").exists():
        env["PATH"] = str(TOOLS_ENV / "bin") + os.pathsep + env.get("PATH", "")
    return env


def _micromamba_platform():
    sysname, mach = platform.system(), platform.machine().lower()
    if sysname == "Linux":
        return "linux-aarch64" if mach in ("aarch64", "arm64") else "linux-64"
    if sysname == "Darwin":
        return "osx-arm64" if mach == "arm64" else "osx-64"
    return None


def get_micromamba(log):
    """Download the static micromamba binary (no privileges required)."""
    plat = _micromamba_platform()
    if plat is None:
        raise RuntimeError("micromamba is available on Linux and macOS only.")
    if not MAMBA_BIN.exists():
        MAMBA_BIN.parent.mkdir(parents=True, exist_ok=True)
        url = f"https://github.com/mamba-org/micromamba-releases/releases/latest/download/micromamba-{plat}"
        log(f"Downloading micromamba ({plat})…")
        tmp = MAMBA_BIN.with_suffix(".part")
        urllib.request.urlretrieve(url, tmp)
        tmp.chmod(0o755)
        tmp.replace(MAMBA_BIN)
    os.environ["MAMBA_ROOT_PREFIX"] = str(APP_HOME / "mamba_root")
    return str(MAMBA_BIN)


LAUNCHER_TAG = "# rnaseqstudio-launcher"
LAUNCH_TOOLS = EXT_TOOLS + ["hisat2-inspect", "fastp", "bgzip", "tabix"]


def make_local_launchers(log):
    """Create small launcher scripts in ~/.local/bin for the installed tools, so they can be used from any
    terminal. Files in ~/.local/bin that were not created by RNAseqStudio are never overwritten."""
    LOCAL_BIN.mkdir(parents=True, exist_ok=True)
    made = []
    for t in dict.fromkeys(LAUNCH_TOOLS):
        target = TOOLS_ENV / "bin" / t
        if not target.exists():
            continue
        dest = LOCAL_BIN / t
        if dest.exists() or dest.is_symlink():
            try:
                if LAUNCHER_TAG not in dest.read_text(errors="ignore")[:200]:
                    log(f"{dest} already exists and was not created by RNAseqStudio: leaving it untouched.")
                    continue
            except (OSError, UnicodeDecodeError):
                continue
        dest.write_text(f'#!/bin/sh\n{LAUNCHER_TAG}\nPATH="{TOOLS_ENV}/bin:$PATH" exec "{target}" "$@"\n')
        dest.chmod(0o755)
        made.append(t)
    return made


def install_bio_tools(log, run_cmd):
    """Download micromamba and create/update the Bioconda environment with the external tools."""
    plat = _micromamba_platform()
    if plat is None:
        raise RuntimeError("Automatic installation of external tools is available on Linux and macOS only. "
                           "On Windows use WSL2; the differential expression and pathway steps work anyway.")
    get_micromamba(log)
    os.environ["MAMBA_ROOT_PREFIX"] = str(APP_HOME / "mamba_root")
    action = "install" if (TOOLS_ENV / "conda-meta").exists() else "create"
    cmd = [str(MAMBA_BIN), action, "-y", "-p", str(TOOLS_ENV), "-c", "conda-forge", "-c", "bioconda"]
    if plat == "osx-arm64":  # some Bioconda packages exist only for osx-64 (Rosetta)
        cmd += ["--platform", "osx-64"]
    log(f"Setting up the Bioconda environment in {TOOLS_ENV} (a few minutes)…")
    run_cmd(cmd + CONDA_PKGS)
    missing = [t for t in EXT_TOOLS if not tool_path(t)]
    if missing:
        raise RuntimeError("Incomplete installation, tools not found: " + ", ".join(missing))
    made = make_local_launchers(log)
    log(f"External tools installed in {TOOLS_ENV}; commands available in {LOCAL_BIN}: " + ", ".join(made))
    if str(LOCAL_BIN) not in os.environ.get("PATH", "").split(os.pathsep):
        log(f"Note: {LOCAL_BIN} is not in the PATH of this session; it will be from the next login "
            f"(or add 'export PATH=\"$HOME/.local/bin:$PATH\"' to ~/.bashrc).")


def cli_install():
    def run_cmd(cmd):
        print("$ " + " ".join(cmd), flush=True)
        subprocess.check_call(cmd, env=tool_env())
    install_bio_tools(lambda m: print(f"[{APP_NAME}] {m}", flush=True), run_cmd)


def _say(msg):
    print(f"[{APP_NAME}] {msg}", flush=True)


def _relaunch():
    """Restart with the same interpreter: required because ~/.local/lib/pythonX.Y/site-packages
    is added to sys.path only if it exists at start-up."""
    args = [sys.executable, os.path.abspath(__file__), *sys.argv[1:]]
    env = dict(os.environ, RNASEQSTUDIO_BOOTSTRAPPED="1")
    if os.name == "nt":
        sys.exit(subprocess.call(args, env=env))
    os.execve(sys.executable, args, env)


def _pip_flags():
    """--user when not inside a virtual environment; --break-system-packages when the distribution
    marks the system Python as 'externally managed' (PEP 668, e.g. Ubuntu 23.04+).
    With --user, packages go to ~/.local only: system files are not modified."""
    flags = []
    if sys.prefix == sys.base_prefix:
        flags.append("--user")
        if (Path(sysconfig.get_paths()["stdlib"]) / "EXTERNALLY-MANAGED").exists():
            flags.append("--break-system-packages")
    return flags


def _ensure_pip():
    if subprocess.call([sys.executable, "-m", "pip", "--version"],
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL) == 0:
        return
    _say("pip not available: installing it in ~/.local with get-pip.py…")
    APP_HOME.mkdir(parents=True, exist_ok=True)
    gp = APP_HOME / "get-pip.py"
    urllib.request.urlretrieve("https://bootstrap.pypa.io/get-pip.py", gp)
    subprocess.check_call([sys.executable, str(gp), "--no-warn-script-location", *_pip_flags()])


def bootstrap():
    """Install missing dependencies into the Python interpreter that launched the program
    (user packages in ~/.local), without virtual environments and without sudo."""
    if os.environ.get("RNASEQSTUDIO_NO_BOOTSTRAP") == "1":
        return
    missing = [m for m in PY_DEPS if importlib.util.find_spec(m) is None]
    if not missing:
        return
    if os.environ.get("RNASEQSTUDIO_BOOTSTRAPPED") == "1":
        sys.exit(f"[{APP_NAME}] Dependencies still missing after installation: {', '.join(missing)}.\n"
                 f"Install them manually: {sys.executable} -m pip install {' '.join(_pip_flags())} "
                 f"{' '.join(PY_DEPS[m] for m in missing)}")
    _say(f"Missing Python dependencies: {', '.join(missing)}. Installing them locally (~/.local)…")
    _ensure_pip()
    subprocess.check_call([sys.executable, "-m", "pip", "install", "--no-warn-script-location", *_pip_flags(),
                           *[PY_DEPS[m] for m in missing]])
    importlib.invalidate_caches()
    _relaunch()


if __name__ == "__main__":
    if "--version" in sys.argv:
        print(f"{APP_NAME} {APP_VERSION}")
        sys.exit(0)
    if "--install-tools" in sys.argv:
        cli_install()
        sys.exit(0)
    bootstrap()

# =============================================================================
# Heavy imports (after bootstrap)
# =============================================================================
import gzip
import json
import re
import shlex
import signal
import inspect
import traceback
from datetime import datetime

import numpy as np
import pandas as pd

from PySide6.QtCore import (Qt, QThread, Signal, QAbstractTableModel, QModelIndex,
                            QSortFilterProxyModel, QSettings, QUrl)
from PySide6.QtGui import QPalette, QColor, QFont, QDesktopServices, QBrush, QShortcut, QKeySequence
from PySide6.QtWidgets import (
    QApplication, QMainWindow, QWidget, QVBoxLayout, QHBoxLayout, QGridLayout, QFormLayout,
    QLabel, QLineEdit, QPushButton, QToolButton, QComboBox, QSpinBox, QDoubleSpinBox, QCheckBox,
    QGroupBox, QTableWidget, QTableWidgetItem, QTableView, QHeaderView, QAbstractItemView,
    QListWidget, QListWidgetItem, QStackedWidget, QSplitter, QPlainTextEdit, QProgressBar,
    QFileDialog, QMessageBox, QInputDialog, QTabWidget, QScrollArea, QFrame, QDialog, QRadioButton, QMenu,
)

import matplotlib
matplotlib.use("QtAgg")
from matplotlib.figure import Figure
from matplotlib.colors import ListedColormap
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qtagg import NavigationToolbar2QT as NavToolbar

UP_COLOR, DOWN_COLOR, NS_COLOR = "#c8364b", "#2b6cb0", "#c3c8d0"
GROUP_COLORS = ["#4b3fa8", "#d08a1e", "#1f8a70", "#c8364b", "#2b6cb0", "#8a4fbf", "#5d6d7e", "#b5651d"]
CONTROL_RE = re.compile(r"(?i)(ctrl|control|controllo|ctr|wt|wild|untreat|sham|veh|mock|basal|baseline|normal)")
SYMBOL_COLS = {"gene_name", "symbol", "gene_symbol", "genesymbol", "gene", "name", "hgnc_symbol",
               "external_gene_name", "mgi_symbol"}
DEFAULT_LIBS = {
    "Human": ["KEGG_2021_Human", "Reactome_2022", "GO_Biological_Process_2023", "MSigDB_Hallmark_2020",
              "WikiPathway_2023_Human", "GO_Molecular_Function_2023", "GO_Cellular_Component_2023"],
    "Mouse": ["KEGG_2019_Mouse", "Reactome_2022", "GO_Biological_Process_2023", "MSigDB_Hallmark_2020",
              "WikiPathways_2019_Mouse", "GO_Molecular_Function_2023", "GO_Cellular_Component_2023"],
}
DEFAULT_CHECKED = {"KEGG_2021_Human", "KEGG_2019_Mouse", "Reactome_2022", "GO_Biological_Process_2023",
                   "MSigDB_Hallmark_2020"}


# =============================================================================
# 1. Background execution
# =============================================================================
class Cancelled(Exception):
    pass


class Runner:
    """Run external commands streaming their output to the log; supports cancellation."""

    def __init__(self, log):
        self.log, self.proc, self.cancelled = log, None, False
        from collections import deque
        self.tail = deque(maxlen=400)

    def run(self, cmd, shell=False, cwd=None):
        if self.cancelled:
            raise Cancelled()
        shown = cmd if isinstance(cmd, str) else " ".join(shlex.quote(str(c)) for c in cmd)
        self.log(f"$ {shown}")
        if shell:
            cmd = ["bash", "-o", "pipefail", "-c", cmd]
        self.proc = subprocess.Popen([str(c) for c in cmd], stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                                     text=True, bufsize=1, env=tool_env(), cwd=cwd, errors="replace",
                                     start_new_session=(os.name != "nt"))
        for line in self.proc.stdout:
            line = line.rstrip()
            if line:
                self.tail.append(line)
                if not line.startswith("UnknownTag="):  # featureCounts noise on corrupted BAM files
                    self.log(line)
        rc = self.proc.wait()
        self.proc = None
        if self.cancelled:
            raise Cancelled()
        if rc != 0:
            raise RuntimeError(f"Command exited with code {rc}: {shown[:200]}")

    def cancel(self):
        self.cancelled = True
        p = self.proc
        if p is not None:
            try:
                if os.name != "nt":
                    os.killpg(p.pid, signal.SIGTERM)
                else:
                    p.terminate()
            except Exception:
                pass


class Task(QThread):
    log = Signal(str)
    progress = Signal(int, int)
    done = Signal(object)
    failed = Signal(str)

    def __init__(self, fn, *args, **kwargs):
        super().__init__()
        self.fn, self.args, self.kwargs = fn, args, kwargs
        self.runner = Runner(self.log.emit)

    def run(self):
        try:
            self.done.emit(self.fn(self, *self.args, **self.kwargs))
        except Cancelled:
            self.failed.emit("Operation cancelled.")
        except Exception as e:
            self.log.emit(traceback.format_exc())
            self.failed.emit(str(e))

    def cancel(self):
        self.runner.cancel()


# =============================================================================
# 2. Utilities: files, samples, counts, annotation
# =============================================================================
def sanitize(s) -> str:
    return re.sub(r"[^A-Za-z0-9._-]+", "_", str(s)).strip("_") or "x"


def _open_text(path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path, encoding="utf-8", errors="replace")


def _sniff_sep(path):
    with _open_text(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            if "\t" in line:
                return "\t"
            return "," if line.count(",") >= line.count(";") else ";"
    return "\t"


def strip_version(ids):
    return pd.Index(ids).astype(str).str.replace(r"^(ENS[A-Z]*\d+)\.\d+$", r"\1", regex=True)


def clean_sample_name(s) -> str:
    n = Path(str(s)).name
    for suf in (".bam", ".sam", ".cram"):
        if n.lower().endswith(suf):
            n = n[: -len(suf)]
    for suf in ("Aligned.sortedByCoord.out", ".sorted", ".sort", ".dedup", ".markdup"):
        if n.endswith(suf):
            n = n[: -len(suf)]
    return n.rstrip("._") or str(s)


def guess_group(name: str) -> str:
    g = re.sub(r"(?i)[_\-. ]*(rep|replicate|r|n)?\d+$", "", name)
    return g or name


FASTQ_EXT = re.compile(r"\.(fastq|fq)(\.gz)?$", re.I)
PAIR_RE = re.compile(r"^(?P<name>.+?)(?:_S\d+)?(?:_L(?P<lane>\d{3}))?[._]R?(?P<read>[12])(?:_001)?\.(?:fastq|fq)(?:\.gz)?$",
                     re.I)


def discover_fastq(folder):
    """Detect FASTQ files (also in sub-folders), pair R1/R2 and group lanes by sample."""
    raw = {}
    for f in sorted(Path(folder).rglob("*")):
        if not f.is_file() or not FASTQ_EXT.search(f.name):
            continue
        m = PAIR_RE.match(f.name)
        if m:
            name, read, lane = m["name"], m["read"], m["lane"] or "000"
        else:
            name = re.sub(r"(_S\d+)?(_L\d{3})?(_001)?$", "", FASTQ_EXT.sub("", f.name))
            read, lane = "1", "000"
        raw.setdefault(name, {"R1": [], "R2": []})["R" + read].append((lane, str(f)))
    out = {}
    for name, d in raw.items():
        r1 = [p for _, p in sorted(d["R1"])]
        r2 = [p for _, p in sorted(d["R2"])]
        if r2 and len(r1) != len(r2):
            r2 = []  # inconsistent pairing: treated as single-end
        if r1:
            out[name] = {"R1": r1, "R2": r2}
    return out


def read_counts(path):
    """Read a count matrix (featureCounts, merged HTSeq/STAR, generic CSV/TSV).
    Returns (counts genes×samples int, lengths or None, symbols or None, non_integer_values)."""
    sep = _sniff_sep(path)
    skip = 0
    with _open_text(path) as fh:
        for line in fh:
            if line.startswith("#"):
                skip += 1
            else:
                break
    df = pd.read_csv(path, sep=sep, skiprows=skip, index_col=0, compression="infer")
    df.index = df.index.astype(str)
    lengths, symbols = None, None
    fc_cols = ["Chr", "Start", "End", "Strand", "Length"]
    if set(fc_cols).issubset(df.columns):
        lengths = pd.to_numeric(df["Length"], errors="coerce")
        df = df.drop(columns=fc_cols)
    num = {}
    for c in df.columns:
        v = pd.to_numeric(df[c], errors="coerce")
        lc = str(c).lower()
        if lc in ("length", "gene_length", "efflength", "effective_length") and v.notna().all():
            lengths = v
        elif v.notna().mean() > 0.9:
            num[c] = v
        elif symbols is None and lc in SYMBOL_COLS:
            symbols = df[c].astype(str)
    counts = pd.DataFrame(num).fillna(0)
    counts = counts[~counts.index.str.startswith(("__", "N_unmapped", "N_multimapping", "N_noFeature", "N_ambiguous"))]
    if counts.shape[1] == 0:
        raise ValueError("No numeric count column found in the file.")
    if (counts.values < 0).any():
        raise ValueError("The file contains negative values: raw counts are required, not normalised/log data.")
    nonint = not np.allclose(counts.values, np.round(counts.values))
    counts = counts.round().astype(np.int64)
    names, seen = [], {}
    for c in counts.columns:
        n = clean_sample_name(c)
        seen[n] = seen.get(n, 0) + 1
        names.append(n if seen[n] == 1 else f"{n}_{seen[n]}")
    counts.columns = names
    if counts.index.duplicated().any():
        counts = counts.groupby(level=0).sum()
    if lengths is not None:
        lengths = lengths[~lengths.index.duplicated()].reindex(counts.index)
    return counts, lengths, symbols, nonint


_ANNOT_CACHE = {}


def parse_gtf_genes(path, log=print):
    log(f"Reading gene annotation from {Path(path).name}…")
    rec = {}
    attr_re = re.compile(r'(\w+) "([^"]*)"')
    with _open_text(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9:
                continue
            a9 = f[8]
            i = a9.find('gene_id "')
            if i < 0:
                continue
            gid = a9[i + 9: a9.find('"', i + 9)]
            if gid in rec and f[2] != "gene":
                continue
            a = dict(attr_re.findall(a9))
            rec[gid] = (a.get("gene_name") or a.get("gene") or a.get("gene_symbol") or gid,
                        a.get("gene_type") or a.get("gene_biotype") or "", f[0])
    return pd.DataFrame.from_dict(rec, orient="index", columns=["gene_name", "gene_type", "chrom"])


def load_annotation(path, log=print):
    path = str(path)
    if path in _ANNOT_CACHE:
        return _ANNOT_CACHE[path]
    if re.search(r"\.gtf(\.gz)?$", path, re.I):
        ann = parse_gtf_genes(path, log)
    else:
        df = pd.read_csv(path, sep=_sniff_sep(path), index_col=0, dtype=str)
        cols = {c.lower(): c for c in df.columns}
        name_col = next((cols[c] for c in ("gene_name", "symbol", "gene_symbol", "external_gene_name",
                                            "hgnc_symbol", "mgi_symbol", "name") if c in cols),
                        df.columns[0] if len(df.columns) else None)
        if name_col is None:
            raise ValueError("The annotation file must have at least one column with the gene symbol.")
        type_col = next((cols[c] for c in ("gene_type", "gene_biotype", "biotype") if c in cols), None)
        ann = pd.DataFrame({"gene_name": df[name_col].values,
                            "gene_type": df[type_col].values if type_col else ""}, index=df.index.astype(str))
    ann = ann[~ann.index.duplicated()]
    _ANNOT_CACHE[path] = ann
    return ann


def map_symbols(ids, ann):
    """Map gene_id -> symbol, handling Ensembl versions (ENSG….12). Returns the ID when unmapped."""
    ids = pd.Index(ids).astype(str)
    fallback = pd.Series(ids, index=ids)
    if ann is None or ann.empty:
        return fallback
    s = pd.Series(ann["gene_name"].values, index=ann.index.astype(str))
    out = s.reindex(ids)
    if out.notna().mean() < 0.5:
        s2 = s.copy()
        s2.index = strip_version(s2.index)
        s2 = s2[~s2.index.duplicated()]
        out2 = pd.Series(s2.reindex(strip_version(ids)).values, index=ids)
        if out2.notna().mean() > out.notna().mean():
            out = out2
    out.index = ids
    return out.fillna(fallback)


def bam_is_paired(bam) -> bool:
    sam = tool_path("samtools") or "samtools"
    p = subprocess.Popen([sam, "view", str(bam)], stdout=subprocess.PIPE, stderr=subprocess.DEVNULL,
                         text=True, env=tool_env())
    paired = False
    try:
        for i, line in enumerate(p.stdout):
            if int(line.split("\t", 2)[1]) & 1:
                paired = True
                break
            if i > 2000:
                break
    finally:
        p.kill()
        p.wait()
    return paired


def bam_contigs(bam):
    out = subprocess.run([tool_path("samtools") or "samtools", "view", "-H", str(bam)], capture_output=True,
                         text=True, env=tool_env())
    return {f[3:] for l in out.stdout.splitlines() if l.startswith("@SQ")
            for f in l.split("\t") if f.startswith("SN:")}


def gtf_contigs(gtf):
    seen = set()
    with _open_text(gtf) as fh:
        for line in fh:
            if line and not line.startswith("#"):
                seen.add(line.split("\t", 1)[0])
    return seen


def _add_chr(c):
    return "chrM" if c == "MT" else (c if c.startswith("chr") else "chr" + c)


def _strip_chr(c):
    return "MT" if c == "chrM" else (c[3:] if c.startswith("chr") else c)


def harmonize_gtf(gtf, bams, outdir, log):
    """Check that GTF and BAM use the same chromosome names; if they differ only by the 'chr'
    prefix (UCSC vs Ensembl/NCBI), write a copy of the GTF with converted names."""
    per_bam = {b: bam_contigs(b) for b in bams}
    empty = [Path(b).name for b, c in per_bam.items() if not c]
    if empty:
        raise RuntimeError("Unreadable header (corrupted BAM or not indexable): " + ", ".join(empty))
    bam_ctg = set.union(*per_bam.values())
    g = gtf_contigs(gtf)
    options = [("no conversion", lambda c: c), ("adding the 'chr' prefix", _add_chr),
               ("removing the 'chr' prefix", _strip_chr)]
    scored = [(len({fn(c) for c in g} & bam_ctg), name, fn) for name, fn in options]
    best = max(scored, key=lambda x: x[0])
    if best[0] == 0:
        raise RuntimeError(
            "No chromosome in common between GTF and BAM: annotation and reference genome do not match.\n"
            f"GTF: {', '.join(sorted(g)[:6])}…\nBAM: {', '.join(sorted(bam_ctg)[:6])}…\n"
            "Use a GTF from the same build and the same source as the genome used for alignment.")
    for b, c in per_bam.items():
        if not ({best[2](x) for x in g} & c):
            log(f"Warning: {Path(b).name} shares no chromosome with the GTF (different reference?).")
    if best[1] == "no conversion":
        return gtf
    log(f"GTF and BAM chromosome names differ: applying {best[1]} to the GTF "
        f"({best[0]} chromosomes in common after conversion).")
    fixed = Path(outdir) / ("annotation_" + re.sub(r"\.gtf(\.gz)?$", "", Path(gtf).name) + "_harmonized.gtf")
    fn = best[2]
    with _open_text(gtf) as fi, open(fixed, "w") as fo:
        for line in fi:
            if line.startswith("#") or "\t" not in line:
                fo.write(line)
            else:
                c, rest = line.split("\t", 1)
                fo.write(fn(c) + "\t" + rest)
    return str(fixed)


def check_bams(bams, log):
    """samtools quickcheck: header and EOF block (detects truncated BAM files)."""
    out = subprocess.run([tool_path("samtools") or "samtools", "quickcheck", "-v", *bams], capture_output=True,
                         text=True, env=tool_env())
    bad = [l.strip() for l in out.stdout.splitlines() if l.strip()]
    missing_bai = [b for b in bams if not (Path(b + ".bai").exists() or Path(b[:-4] + ".bai").exists())]
    if missing_bai:
        log("Note: missing .bai index for " + ", ".join(Path(b).name for b in missing_bai) + " (not required).")
    return bad


def featurecounts_version():
    try:
        out = subprocess.run([tool_path("featureCounts") or "featureCounts", "-v"], capture_output=True,
                             text=True, env=tool_env())
        m = re.search(r"v(\d+)\.(\d+)\.(\d+)", out.stdout + out.stderr)
        return tuple(int(x) for x in m.groups()) if m else (0, 0, 0)
    except Exception:
        return (0, 0, 0)


def nospace(path):
    """Return an equivalent path without spaces (HISAT2 and STAR do not handle spaces):
    symbolic link to the folder, or to the file when the file name itself contains spaces."""
    import hashlib
    p = str(path)
    if " " not in p:
        return p
    real = Path(p)
    links = APP_HOME / "links"
    links.mkdir(parents=True, exist_ok=True)
    if " " in real.name:
        target = real.resolve()
        link = links / (hashlib.sha1(str(target).encode()).hexdigest()[:10] + "_" + sanitize(real.name))
        if not link.is_symlink():
            link.symlink_to(target)
        return str(link)
    parent = real.parent.resolve()
    link = links / hashlib.sha1(str(parent).encode()).hexdigest()[:12]
    if not link.is_symlink():
        link.symlink_to(parent, target_is_directory=True)
    return str(link / real.name)


def ensure_uncompressed(path, dest_dir, log):
    """HISAT2/STAR cannot read compressed FASTA: decompress once into the index folder."""
    path = Path(path)
    if not path.name.endswith(".gz"):
        return str(path)
    out = Path(dest_dir) / path.name[:-3]
    if not out.exists():
        log(f"Decompressing {path.name}…")
        with gzip.open(path, "rb") as fi, open(out, "wb") as fo:
            shutil.copyfileobj(fi, fo, 16 << 20)
    return str(out)


def bh(p, m=None):
    p = np.asarray(p, float)
    n = len(p)
    if n == 0:
        return p
    m = max(m or n, n)
    order = np.argsort(p)
    ranked = p[order] * m / np.arange(1, n + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty(n)
    out[order] = np.minimum(ranked, 1.0)
    return out


def read_gmt(path):
    gs = {}
    with _open_text(path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) >= 3:
                gs[f[0]] = [g for g in f[2:] if g]
    return gs


def write_gmt(gs, path):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as fh:
        for k, v in gs.items():
            fh.write("\t".join([k, ""] + list(v)) + "\n")


def load_library(name, organism, log):
    cache = GENESET_CACHE / f"{organism}__{sanitize(name)}.gmt"
    if cache.exists():
        return read_gmt(cache)
    import gseapy as gp
    log(f"Downloading library {name} ({organism}) from Enrichr…")
    gs = gp.get_library(name=name, organism=organism)
    write_gmt(gs, cache)
    return gs


# =============================================================================
# 3. Analysis steps (run in worker threads)
# =============================================================================
def task_install_tools(task):
    install_bio_tools(task.log.emit, task.runner.run)
    return True


def task_align(task, samples, cfg):
    log, R = task.log.emit, task.runner
    out, T, aligner = Path(cfg["outdir"]), int(cfg["threads"]), cfg["aligner"]
    need = ["samtools", "hisat2" if aligner == "HISAT2" else "STAR"] + (["fastp"] if cfg["trim"] else [])
    for t in need:
        if not tool_path(t):
            raise RuntimeError(f"'{t}' not found. Install it from the Environment page.")
    idx = cfg["index"].strip()
    if aligner == "HISAT2":
        if idx:
            idx = re.sub(r"\.\d+\.ht2l?$", "", idx)
        else:
            if not cfg["fasta"]:
                raise ValueError("Provide the genome FASTA or an existing HISAT2 index.")
            idxdir = out / "00_index" / "hisat2"
            idxdir.mkdir(parents=True, exist_ok=True)
            idx = str(idxdir / "genome")
            if not (Path(idx + ".1.ht2").exists() or Path(idx + ".1.ht2l").exists()):
                fa = ensure_uncompressed(cfg["fasta"], idxdir, log)
                log("Building the HISAT2 index (human genome: about 1 hour and ~8 GB of RAM)…")
                R.run(["hisat2-build", "-p", T, nospace(fa), nospace(idx)])
    else:
        if not idx:
            if not cfg["fasta"]:
                raise ValueError("Provide the genome FASTA or an existing STAR index.")
            idxdir = out / "00_index" / "star"
            idxdir.mkdir(parents=True, exist_ok=True)
            idx = str(idxdir)
            if not (idxdir / "SAindex").exists():
                fa = ensure_uncompressed(cfg["fasta"], idxdir, log)
                log("Building the STAR index (human genome: ~32 GB of RAM required)…")
                sidx = nospace(idx)
                cmd = ["STAR", "--runMode", "genomeGenerate", "--runThreadN", T, "--genomeDir", sidx,
                       "--genomeFastaFiles", nospace(fa), "--outFileNamePrefix", sidx + "/"]
                if cfg["gtf"]:
                    gtf = ensure_uncompressed(cfg["gtf"], idxdir, log)
                    cmd += ["--sjdbGTFfile", nospace(gtf), "--sjdbOverhang", cfg["overhang"]]
                R.run(cmd)

    bamdir = out / "02_alignment"
    bamdir.mkdir(parents=True, exist_ok=True)
    trimdir = out / "01_trimming"
    q = shlex.quote
    extra = cfg["extra"].strip()
    sidx = nospace(idx)

    def trim(name, r1, r2, paired):
        trimdir.mkdir(parents=True, exist_ok=True)
        t1, t2 = [], []
        for j, a in enumerate(r1):
            o1 = trimdir / f"{name}_L{j + 1}_R1.trim.fastq.gz"
            o2 = trimdir / f"{name}_L{j + 1}_R2.trim.fastq.gz"
            js = trimdir / f"{name}_L{j + 1}.fastp.json"  # written by fastp only when finished
            if js.exists() and o1.exists() and (not paired or o2.exists()):
                log(f"Trimming of {name} (lane {j + 1}) already present: reusing it.")
            else:
                js.unlink(missing_ok=True)
                cmd = ["fastp", "-w", min(T, 16), "-i", nospace(a), "-o", nospace(o1),
                       "-h", nospace(trimdir / f"{name}_L{j + 1}.fastp.html"), "-j", nospace(js)]
                if paired:
                    cmd += ["-I", nospace(r2[j]), "-O", nospace(o2), "--detect_adapter_for_pe"]
                try:
                    R.run(cmd)
                except Exception:
                    for f in (o1, o2, js):
                        Path(f).unlink(missing_ok=True)
                    raise
            t1.append(str(o1))
            if paired:
                t2.append(str(o2))
        return t1, t2

    def align(name, r1, r2, paired, bam):
        part = bamdir / f"{name}.part.bam"  # renamed to .bam only when complete
        for f in bamdir.glob(f"{name}.part.bam*"):
            f.unlink()
        r1, r2 = [nospace(x) for x in r1], [nospace(x) for x in r2]
        spart = nospace(part)
        if aligner == "HISAT2":
            summ = bamdir / f"{name}.hisat2.summary.txt"
            reads = (f"-1 {q(','.join(r1))} -2 {q(','.join(r2))}" if paired else f"-U {q(','.join(r1))}")
            cmd = (f"hisat2 -p {T} -x {q(sidx)} {reads} --new-summary --summary-file {q(nospace(summ))} {extra} "
                   f"| samtools sort -@ {max(1, T // 2)} -m 1G -T {q(spart + '.tmp')} -o {q(spart)} - "
                   f"&& samtools index {q(spart)}")
            R.run(cmd, shell=True)
            if summ.exists():
                for line in summ.read_text().splitlines():
                    if "overall" in line.lower():
                        log(f"[{name}] {line.strip()}")
        else:
            prefix = nospace(bamdir / f"{name}.")
            cmd = ["STAR", "--runThreadN", T, "--genomeDir", sidx, "--readFilesIn", ",".join(r1)]
            if paired:
                cmd.append(",".join(r2))
            if r1[0].endswith(".gz"):
                cmd += ["--readFilesCommand", "gzip", "-cd"]
            cmd += ["--outSAMtype", "BAM", "SortedByCoordinate", "--outSAMattributes", "NH", "HI", "AS", "nM",
                    "--outFileNamePrefix", prefix] + shlex.split(extra)
            R.run(cmd)
            os.replace(prefix + "Aligned.sortedByCoord.out.bam", part)
            R.run(["samtools", "index", spart])
            lf = Path(prefix + "Log.final.out")
            if lf.exists():
                for line in lf.read_text().splitlines():
                    if "Uniquely mapped reads %" in line or "% of reads mapped to multiple loci" in line:
                        log(f"[{name}] {line.strip()}")
        os.replace(str(part) + ".bai", str(bam) + ".bai")
        os.replace(part, bam)

    bams, failed, any_paired, n = [], {}, False, len(samples)
    for i, (name, fq) in enumerate(samples.items()):
        task.progress.emit(i, n)
        log(f"━━ Sample {i + 1}/{n}: {name}")
        bam = bamdir / f"{name}.bam"
        r1, r2 = list(fq["R1"]), list(fq["R2"])
        paired = bool(r2)
        if bam.exists() and Path(str(bam) + ".bai").exists():
            log(f"{bam.name} already complete: skipping (delete it to redo the alignment).")
            bams.append(str(bam))
            any_paired |= paired
            continue
        try:
            if cfg["trim"]:
                r1, r2 = trim(name, r1, r2, paired)
            align(name, r1, r2, paired, bam)
            bams.append(str(bam))
            any_paired |= paired
        except Cancelled:
            raise
        except Exception as e:
            for f in bamdir.glob(f"{name}.part.bam*"):
                f.unlink(missing_ok=True)
            if not cfg.get("skip_errors", True):
                raise
            reason = next((l.strip() for l in reversed(R.tail) if re.search(r"(?i)error|truncated|abort", l)), str(e))
            failed[name] = reason
            log(f"✖ {name} not aligned ({reason}). Continuing with the next samples.")
    task.progress.emit(n, n)
    if failed:
        pd.DataFrame({"sample": list(failed), "error": list(failed.values())}).to_csv(
            bamdir / "failed_samples.tsv", sep="\t", index=False)
        log(f"Samples not aligned ({len(failed)}): {', '.join(failed)} — details in "
            f"{bamdir / 'failed_samples.tsv'}")
    if not bams:
        raise RuntimeError("No sample was aligned successfully.")
    return {"bams": bams, "paired": any_paired, "gtf": cfg["gtf"], "failed": failed}


def task_featurecounts(task, bams, cfg):
    log, R = task.log.emit, task.runner
    for t in ("featureCounts", "samtools"):
        if not tool_path(t):
            raise RuntimeError(f"'{t}' not found. Install it from the Environment page.")
    if not cfg["gtf"]:
        raise ValueError("Provide the GTF annotation file.")
    out = Path(cfg["outdir"]) / "03_counts"
    out.mkdir(parents=True, exist_ok=True)
    skip_bad = cfg.get("skip_bad", True)
    excluded = {}
    log("Checking BAM integrity (samtools quickcheck)…")
    bad = check_bams(bams, log)
    if bad:
        if not skip_bad:
            raise RuntimeError("Corrupted or truncated BAM files:\n  " + "\n  ".join(bad) +
                               "\nRegenerate (or copy again) or remove them from the list.")
        for b in bad:
            excluded[b] = "truncated or invalid header (quickcheck)"
            log(f"✖ Excluding {Path(b).name}: {excluded[b]}")
        bams = [b for b in bams if b not in bad]
        if not bams:
            raise RuntimeError("All BAM files are damaged.")
    gtf = harmonize_gtf(cfg["gtf"], bams, out, log)
    paired = cfg["paired"]
    if paired is None:
        paired = bam_is_paired(bams[0])
        log(f"Library type detected: {'paired-end' if paired else 'single-end'}")
    raw = out / "featureCounts.txt"
    sig_file = out / "featureCounts.params.json"

    def signature(bl):
        return {"gtf": str(gtf), "strand": cfg["strand"], "feature": cfg["feature"], "attribute": cfg["attribute"],
                "mapq": cfg["mapq"], "paired": bool(paired),
                "bams": [[str(b), Path(b).stat().st_size, int(Path(b).stat().st_mtime)] for b in bl]}

    reuse = False
    if raw.exists() and Path(str(raw) + ".summary").exists() and sig_file.exists():
        try:
            reuse = json.loads(sig_file.read_text()) == signature(bams)
        except Exception:
            reuse = False
    if reuse:
        log("Counts already computed with the same BAM files and parameters: reusing featureCounts.txt "
            "(delete it to recompute).")
    while not reuse:
        for f in (raw, Path(str(raw) + ".summary"), sig_file):
            f.unlink(missing_ok=True)
        by_name = {Path(nospace(b)).name: b for b in bams}
        cmd = ["featureCounts", "-T", cfg["threads"], "-a", nospace(gtf), "-o", nospace(raw), "-s", cfg["strand"],
               "-t", cfg["feature"], "-g", cfg["attribute"], "-Q", cfg["mapq"]]
        if paired:
            cmd += ["-p", "-B", "-C"]
            if featurecounts_version() >= (2, 0, 2):
                cmd.append("--countReadPairs")
        R.tail.clear()
        R.run(cmd + [nospace(b) for b in bams])
        tail = list(R.tail)
        if raw.exists() and not any("GZIP ERROR" in l or l.startswith("ERROR") for l in tail):
            sig_file.write_text(json.dumps(signature(bams)))
            break
        m = next((re.search(r"Process BAM file (\S+?)\.\.\.", l) for l in reversed(tail)
                  if "Process BAM file" in l), None)
        current = by_name.get(m.group(1)) if m else None
        err = next((l.strip() for l in reversed(tail) if "ERROR" in l), "no output produced")
        if current and skip_bad and len(bams) > 1:
            excluded[current] = f"damaged internally ({err})"
            log(f"✖ featureCounts stopped on {Path(current).name} ({err}): excluding it and restarting "
                f"with the remaining {len(bams) - 1} BAM.")
            bams = [b for b in bams if b != current]
            continue
        msg = f"featureCounts stopped ({err})"
        if current:
            msg += (f" while reading {Path(current).name}.\nThe BAM is damaged internally: "
                    f"check it with 'samtools view -c {current} > /dev/null', then regenerate or exclude it.")
        raise RuntimeError(msg)
    if excluded:
        pd.DataFrame({"bam": list(excluded), "reason": list(excluded.values())}).to_csv(
            out / "excluded_bams.tsv", sep="\t", index=False)

    counts, lengths, _, _ = read_counts(raw)
    counts_path = out / "counts_matrix.tsv"
    counts.to_csv(counts_path, sep="\t", index_label="gene_id")
    if lengths is not None:
        rpk = counts.div(lengths.replace(0, np.nan) / 1e3, axis=0).fillna(0)
        (rpk.div(rpk.sum(axis=0), axis=1) * 1e6).round(3).to_csv(out / "tpm_matrix.tsv", sep="\t",
                                                                  index_label="gene_id")
    ann = load_annotation(gtf, log)
    ann_path = out / "gene_annotation.tsv"
    ann.to_csv(ann_path, sep="\t", index_label="gene_id")

    summary = None
    sfile = Path(str(raw) + ".summary")
    if sfile.exists():
        s = pd.read_csv(sfile, sep="\t", index_col=0)
        s.columns = [clean_sample_name(c) for c in s.columns]
        s = s.loc[:, s.sum() > 0]
        tot = s.sum()
        summary = pd.DataFrame({"Sample": list(s.columns), "Total": tot.values,
                                "Assigned reads": s.loc["Assigned"].values,
                                "% assigned": (100 * s.loc["Assigned"] / tot).round(1).values})
        for row in ("Unassigned_NoFeatures", "Unassigned_Ambiguity", "Unassigned_MultiMapping"):
            if row in s.index:
                summary[row.replace("Unassigned_", "% ")] = (100 * s.loc[row] / tot).round(1).values
        summary.to_csv(out / "assignment_summary.tsv", sep="\t", index=False)
        zero = summary[summary["Assigned reads"] == 0]
        if len(zero) == len(summary):
            raise RuntimeError("No read assigned in any sample: check that GTF and BAM come "
                               "from the same genome and that -t/-g match the GTF (exon / gene_id).")
        low = summary[summary["% assigned"] < 40]
        if len(low):
            log("Warning: assignment < 40% in " + ", ".join(low["Sample"]) +
                ". Check the strandedness (-s) or the annotation.")
    log(f"Count matrix: {counts.shape[0]} genes × {counts.shape[1]} samples → {counts_path}")
    return {"counts_path": str(counts_path), "annotation_path": str(ann_path), "summary": summary,
            "excluded": excluded}


def _construct(cls, **kw):
    """Instantiate PyDESeq2 classes passing only the arguments supported by the installed version."""
    params = inspect.signature(cls).parameters
    if any(p.kind == p.VAR_KEYWORD for p in params.values()):
        return cls(**kw)
    return cls(**{k: v for k, v in kw.items() if k in params})


def _mor_normalize(X):
    """Median-of-ratios normalisation (DESeq2) on a samples × genes matrix."""
    lg = np.log(X.where(X > 0))
    geo = lg.mean(axis=0)
    ok = lg.notna().all(axis=0)
    if ok.sum() == 0:
        return X.astype(float)
    sf = np.exp((lg.loc[:, ok] - geo[ok]).median(axis=1))
    return X.div(sf, axis=0)


def _blind_fit(DeseqDataSet, kw, log):
    """Without replicates: dispersions estimated with design ~1 (all samples treated as replicates),
    then fit of the full model with those dispersions fixed."""
    if "design" not in kw:
        raise RuntimeError("requires PyDESeq2 ≥ 0.5")
    log("Blind dispersion estimation (design ~1)…")
    kb = dict(kw, design="~1", refit_cooks=False)
    b = _construct(DeseqDataSet, **kb)
    b.fit_size_factors()
    b.fit_genewise_dispersions()
    b.fit_dispersion_trend()
    b.fit_dispersion_prior()
    b.fit_MAP_dispersions()
    d = _construct(DeseqDataSet, **dict(kw, refit_cooks=False))
    d.fit_size_factors()
    for col in b.var.columns:
        d.var[col] = b.var[col].values
    for attr in ("non_zero_idx", "non_zero_genes"):
        if hasattr(b, attr):
            setattr(d, attr, getattr(b, attr))
    for k, v in b.uns.items():
        d.uns.setdefault(k, v)
    log("GLM fit with fixed dispersions…")
    d.fit_LFC()
    return d


def task_deseq(task, counts, groups, batches, ref, case, cfg):
    log = task.log.emit
    from pydeseq2.dds import DeseqDataSet
    from pydeseq2.ds import DeseqStats

    samples = [s for s in counts.columns if s in groups]
    if cfg["only_contrast"]:
        samples = [s for s in samples if groups[s] in (ref, case)]
    n_ref = sum(groups[s] == ref for s in samples)
    n_case = sum(groups[s] == case for s in samples)
    warnings = []
    if n_ref < 2 or n_case < 2:
        warnings.append(f"Insufficient replicates (control '{ref}': {n_ref}, case '{case}': {n_case}). "
                        "Biological variability cannot be estimated for this comparison: results should be considered "
                        "exploratory and require experimental validation.")
        log("⚠ " + warnings[-1])
    levels = [ref, case] + sorted(set(groups[s] for s in samples) - {ref, case})
    tok = {g: f"L{i:02d}{re.sub(r'[^A-Za-z0-9]', '', g)}" for i, g in enumerate(levels)}
    meta = pd.DataFrame({"condition": [tok[groups[s]] for s in samples]}, index=samples)
    use_batch = cfg["use_batch"] and all(batches.get(s) for s in samples) and len({batches[s] for s in samples}) > 1
    if cfg["use_batch"] and not use_batch:
        log("Batch covariate ignored: missing values or a single level.")
    if use_batch:
        meta["batch"] = ["B" + re.sub(r"[^A-Za-z0-9]", "", str(batches[s])) for s in samples]
        D = pd.get_dummies(meta[["batch", "condition"]], drop_first=True).astype(float)
        D.insert(0, "intercept", 1.0)
        if np.linalg.matrix_rank(D.values) < D.shape[1]:
            raise ValueError("Batch and group are confounded (design not of full rank): "
                             "remove the batch covariate or revise the assignment.")

    c = counts[samples]
    min_n = cfg["min_samples"] or min(n_ref, n_case)
    keep = (c >= cfg["min_count"]).sum(axis=1) >= min_n
    c = c[keep]
    log(f"Filter: {int(keep.sum())}/{len(keep)} genes with ≥{cfg['min_count']} reads in ≥{min_n} samples.")
    log(f"Model: {'~ batch + condition' if use_batch else '~ condition'} — "
        f"{len(samples)} samples, contrast '{case}' vs '{ref}'.")
    X = c.T.astype(np.int64)

    Dm = pd.get_dummies(meta[[c for c in ("batch", "condition") if c in meta.columns]], drop_first=True).astype(float)
    Dm.insert(0, "intercept", 1.0)
    resid_df = len(samples) - int(np.linalg.matrix_rank(Dm.values))
    mode = "standard" if resid_df > 0 else "blind"
    if mode == "blind":
        warnings.append("No replicates in the model (0 residual degrees of freedom): dispersion estimated 'blindly', "
                        "treating all samples as replicates (as in DESeq2 ≤ 1.20). Very conservative p-values.")
        log("⚠ " + warnings[-1])
    elif min(n_ref, n_case) < 2:
        log(f"Dispersion is estimated from the other replicated groups ({resid_df} residual degrees of freedom).")

    params = inspect.signature(DeseqDataSet).parameters
    kw = dict(counts=X, metadata=meta, refit_cooks=True, quiet=True, n_cpus=int(cfg["threads"]))
    inference = None
    try:
        from pydeseq2.default_inference import DefaultInference
        inference = DefaultInference(n_cpus=int(cfg["threads"]))
        kw["inference"] = inference
    except Exception:
        pass
    if "design" in params:
        kw["design"] = "~batch + condition" if use_batch else "~condition"
    else:
        kw["design_factors"] = ["batch", "condition"] if use_batch else "condition"
        kw["ref_level"] = ["condition", tok[ref]]
    dds, ds = None, None
    if mode == "standard":
        dds = _construct(DeseqDataSet, **kw)
        log("Estimating size factors, dispersions and GLM fit (PyDESeq2)…")
        dds.deseq2()
    else:
        try:
            dds = _blind_fit(DeseqDataSet, kw, log)
        except Exception as e:
            mode = "exploratory"
            warnings.append(f"Blind estimation failed ({e}): only fold changes on normalised counts were computed, "
                            "without p-values.")
            log("⚠ " + warnings[-1])

    lfc_col = "log2FoldChange"
    if mode == "exploratory":
        norm = _mor_normalize(X).T
        mref = norm[[s for s in samples if groups[s] == ref]].mean(axis=1)
        mcase = norm[[s for s in samples if groups[s] == case]].mean(axis=1)
        res = pd.DataFrame({"baseMean": norm.mean(axis=1), "log2FoldChange": np.log2((mcase + 0.5) / (mref + 0.5)),
                            "lfcSE": np.nan, "stat": np.nan, "pvalue": np.nan, "padj": np.nan})
    else:
        skw = dict(dds=dds, contrast=["condition", tok[case], tok[ref]], alpha=cfg["alpha"],
                   cooks_filter=(mode == "standard"), independent_filter=True, quiet=True)
        if inference is not None:
            skw["inference"] = inference
        ds = _construct(DeseqStats, **skw)
        ds.summary()
        res = ds.results_df.copy()

    if cfg["shrink"] and ds is not None:
        try:
            coeffs = list(dds.varm["LFC"].columns)
            coeff = next(cn for cn in coeffs if tok[case] in cn and "condition" in cn)
            ds.lfc_shrink(coeff=coeff)
            res["log2FoldChange_shrunk"] = ds.results_df["log2FoldChange"].values
            lfc_col = "log2FoldChange_shrunk"
            log("LFC shrinkage (apeGLM) applied.")
        except Exception as e:
            log(f"Shrinkage not applied: {e}")

    if dds is not None:
        norm = pd.DataFrame(np.asarray(dds.layers["normed_counts"]).T, index=X.columns, columns=X.index)
    try:
        if dds is None:
            raise RuntimeError("exploratory mode")
        if hasattr(dds, "vst_fit"):
            dds.vst_fit(use_design=False)
            v = dds.vst_transform()
        else:
            dds.vst(use_design=False)
            v = dds.layers["vst_counts"]
        vst = pd.DataFrame(np.asarray(v).T, index=X.columns, columns=X.index)
    except Exception as e:
        log(f"VST not available ({e}): using log2(normalised counts + 1).")
        vst = np.log2(norm + 1)

    ann = None
    if cfg.get("annotation"):
        try:
            ann = load_annotation(cfg["annotation"], log)
        except Exception as e:
            log(f"Annotation not loaded: {e}")
    symbols = map_symbols(res.index, ann)
    if cfg.get("symbols") is not None and (ann is None):
        symbols = cfg["symbols"].reindex(res.index).fillna(symbols)
    res.insert(0, "symbol", symbols.values)
    if ann is not None and "gene_type" in ann.columns:
        gt = ann["gene_type"].reindex(res.index)
        if gt.isna().mean() > 0.5:
            a2 = ann["gene_type"].copy()
            a2.index = strip_version(a2.index)
            a2 = a2[~a2.index.duplicated()]
            gt = pd.Series(a2.reindex(strip_version(res.index)).values, index=res.index)
        res.insert(1, "gene_type", gt.fillna("").values)
    ref_s = [s for s in samples if groups[s] == ref]
    case_s = [s for s in samples if groups[s] == case]
    res[f"mean_norm_{ref}"] = norm[ref_s].mean(axis=1).reindex(res.index).values
    res[f"mean_norm_{case}"] = norm[case_s].mean(axis=1).reindex(res.index).values
    if mode == "exploratory":
        sig = (res[lfc_col].abs() >= cfg["lfc"]) & (res["baseMean"] >= cfg["min_count"])
    else:
        sig = (res["padj"] < cfg["alpha"]) & (res[lfc_col].abs() >= cfg["lfc"])
    res["direction"] = np.where(sig & (res[lfc_col] > 0), "UP", np.where(sig & (res[lfc_col] < 0), "DOWN", ""))
    res = res.sort_values("padj", na_position="last")
    res.index.name = "gene_id"

    outdir = Path(cfg["outdir"]) / "04_DE" / f"{sanitize(case)}_vs_{sanitize(ref)}"
    outdir.mkdir(parents=True, exist_ok=True)
    full = res.join(norm.add_prefix("norm_"))
    full.to_csv(outdir / "DE_results_full.tsv", sep="\t")
    res.to_csv(outdir / "DE_results.tsv", sep="\t")
    res[res["direction"] != ""].to_csv(outdir / "DE_significant.tsv", sep="\t")
    norm.to_csv(outdir / "normalized_counts.tsv", sep="\t", index_label="gene_id")
    vst.to_csv(outdir / "vst_counts.tsv", sep="\t", index_label="gene_id")
    md = pd.DataFrame({"group": [groups[s] for s in samples]}, index=pd.Index(samples, name="sample"))
    if use_batch:
        md["batch"] = [batches[s] for s in samples]
    md.to_csv(outdir / "sample_metadata.tsv", sep="\t")
    params_out = {k: v for k, v in cfg.items() if k not in ("symbols",)}
    params_out.update(estimation_mode=mode, warnings=warnings, control=ref, case=case, samples=samples, genes_tested=int(len(res)),
                      pydeseq2=_pkg_version("pydeseq2"), date=datetime.now().isoformat(timespec="seconds"))
    (outdir / "parameters.json").write_text(json.dumps(params_out, indent=2, ensure_ascii=False, default=str))
    n_up, n_down = int((res["direction"] == "UP").sum()), int((res["direction"] == "DOWN").sum())
    log(f"Results: {n_up} genes UP and {n_down} DOWN (padj < {cfg['alpha']}, |log2FC| ≥ {cfg['lfc']}) → {outdir}")
    return {"res": res, "vst": vst, "groups": {s: groups[s] for s in samples},
            "batches": {s: batches.get(s, "") for s in samples} if use_batch else {},
            "ref": ref, "case": case, "lfc_col": lfc_col, "outdir": str(outdir), "cfg": cfg,
            "n_up": n_up, "n_down": n_down, "warnings": warnings, "mode": mode}


def ora(genes, background, gene_sets, min_size, max_size):
    """Hypergeometric over-representation analysis with universe = tested genes."""
    from scipy.stats import hypergeom
    bg = set(background)
    genes = set(genes) & bg
    N, n = len(bg), len(genes)
    rows, tested = [], 0
    for term, members in gene_sets.items():
        m = set(members) & bg
        K = len(m)
        if K < min_size or K > max_size:
            continue
        tested += 1
        hit = genes & m
        k = len(hit)
        if k == 0:
            continue
        p = hypergeom.sf(k - 1, N, K, n)
        exp = n * K / N if N else np.nan
        rows.append((term, f"{k}/{K}", k, K, p, k / exp if exp else np.nan, ";".join(sorted(hit))))
    df = pd.DataFrame(rows, columns=["Term", "Overlap", "k", "Size", "P-value", "Fold enrichment", "Genes"])
    if not df.empty:
        df.insert(5, "FDR", bh(df["P-value"].values, tested))
    return df


def task_pathway(task, res, cfg):
    log = task.log.emit
    T = int(cfg["threads"])
    df = res.copy()
    sym = df["symbol"].astype(str) if "symbol" in df.columns else pd.Series(strip_version(df.index), index=df.index)
    df["_g"] = sym.str.strip().str.upper().values
    df = df[(df["_g"] != "") & (df["_g"] != "NAN")]
    lfc = cfg["lfc_col"] if cfg["lfc_col"] in df.columns else "log2FoldChange"

    sets = {}
    for lib in cfg["libraries"]:
        try:
            sets[lib] = load_library(lib, cfg["organism"], log)
        except Exception as e:
            log(f"Library {lib} not available: {e}")
    for g in cfg["gmts"]:
        sets[Path(g).name.replace(".gmt", "")] = read_gmt(g)
    if not sets:
        raise ValueError("No gene-set library available (check the connection or add GMT files).")
    sets = {lib: {t: [x.upper() for x in m] for t, m in gs.items()} for lib, gs in sets.items()}

    outdir = Path(cfg["outdir"])
    outdir.mkdir(parents=True, exist_ok=True)
    bg = df["_g"][df["baseMean"] > 0].unique() if "baseMean" in df.columns else df["_g"].unique()
    ora_df, gsea_df, gsea_obj, ranking = pd.DataFrame(), pd.DataFrame(), {}, None

    if cfg["ora"]:
        sig = (df["padj"] < cfg["alpha"]) & (df[lfc].abs() >= cfg["lfc"])
        up = df.loc[sig & (df[lfc] > 0), "_g"].unique()
        down = df.loc[sig & (df[lfc] < 0), "_g"].unique()
        log(f"ORA: {len(up)} UP genes, {len(down)} DOWN; universe {len(bg)} tested genes.")
        frames = []
        for lib, gs in sets.items():
            for direction, genes in (("UP", up), ("DOWN", down)):
                if len(genes) < 3:
                    continue
                r = ora(genes, bg, gs, cfg["min_size"], cfg["max_size"])
                if not r.empty:
                    r.insert(0, "Library", lib)
                    r.insert(0, "Direction", direction)
                    frames.append(r)
        if frames:
            ora_df = pd.concat(frames, ignore_index=True).sort_values("FDR")
            ora_df.to_csv(outdir / "ORA_results.tsv", sep="\t", index=False)
        log(f"ORA: {int((ora_df['FDR'] < 0.05).sum()) if len(ora_df) else 0} terms with FDR < 0.05.")

    if cfg["gsea"]:
        import gseapy as gp
        if cfg["metric"] == "stat" and "stat" in df.columns:
            score = df["stat"]
        elif cfg["metric"] == "lfc":
            score = df[lfc]
        else:
            score = np.sign(df[lfc]) * -np.log10(df["pvalue"].clip(lower=1e-300))
        r = pd.DataFrame({"g": df["_g"].values, "s": pd.to_numeric(score, errors="coerce").values})
        r = r.replace([np.inf, -np.inf], np.nan).dropna()
        r = r.reindex(r["s"].abs().sort_values(ascending=False).index).drop_duplicates("g")
        ranking = r.set_index("g")["s"].sort_values(ascending=False)
        ranking.to_csv(outdir / "GSEA_ranking.rnk", sep="\t", header=False)
        frames = []
        for lib, gs in sets.items():
            log(f"Preranked GSEA on {lib} ({cfg['perm']} permutations)…")
            kw = dict(rnk=ranking, gene_sets=gs, threads=T, permutation_num=cfg["perm"], min_size=cfg["min_size"],
                      max_size=cfg["max_size"], seed=42, outdir=None, no_plot=True, verbose=False)
            try:
                pre = gp.prerank(**kw)
            except TypeError:
                kw.pop("no_plot")
                pre = gp.prerank(**kw)
            d = pre.res2d.copy()
            if "Term" not in d.columns:
                d = d.reset_index().rename(columns={d.index.name or "index": "Term"})

            def col(*names):
                return next((d[n] for n in names if n in d.columns), pd.Series([np.nan] * len(d)))
            g = pd.DataFrame({
                "Library": lib, "Term": col("Term").values,
                "NES": pd.to_numeric(col("NES", "nes"), errors="coerce").values,
                "ES": pd.to_numeric(col("ES", "es"), errors="coerce").values,
                "NOM p": pd.to_numeric(col("NOM p-val", "pval"), errors="coerce").values,
                "FDR": pd.to_numeric(col("FDR q-val", "fdr"), errors="coerce").values,
                "Tag %": col("Tag %", "tag %").values,
                "Leading edge": col("Lead_genes", "ledge_genes", "lead_genes").values})
            frames.append(g)
            gsea_obj[lib] = getattr(pre, "results", {})
        if frames:
            gsea_df = pd.concat(frames, ignore_index=True).sort_values(["FDR", "NOM p"])
            gsea_df.to_csv(outdir / "GSEA_results.tsv", sep="\t", index=False)
        log(f"GSEA: {int((gsea_df['FDR'] < 0.25).sum()) if len(gsea_df) else 0} terms with FDR < 0.25.")
    return {"ora": ora_df, "gsea": gsea_df, "gsea_obj": gsea_obj, "ranking": ranking, "outdir": str(outdir)}


def _pkg_version(name):
    try:
        from importlib.metadata import version
        return version(name)
    except Exception:
        return None


# =============================================================================
# 4. Plots
# =============================================================================
def _ordered_samples(groups, ref=None, case=None):
    order = [g for g in (ref, case) if g] + sorted(set(groups.values()) - {ref, case})
    return sorted(groups, key=lambda s: (order.index(groups[s]) if groups[s] in order else 99, s))


def _group_colors(groups, ref=None, case=None):
    order = [g for g in (ref, case) if g] + sorted(set(groups.values()) - {ref, case})
    return {g: GROUP_COLORS[i % len(GROUP_COLORS)] for i, g in enumerate(order)}


def _style(ax):
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    ax.tick_params(labelsize=9)


def plot_pca(fig, vst, groups, batches, ref, case, ntop=500):
    ax = fig.add_subplot(111)
    v = vst[[s for s in groups if s in vst.columns]]
    v = v.loc[v.var(axis=1).sort_values(ascending=False).index[:ntop]]
    X = (v.T - v.T.mean()).values
    U, S, _ = np.linalg.svd(X, full_matrices=False)
    pcs, var = U * S, S ** 2 / np.sum(S ** 2)
    colors = _group_colors(groups, ref, case)
    markers = ["o", "s", "^", "D", "v", "P", "X", "*"]
    blevels = sorted({batches.get(s, "") for s in v.columns}) if batches else [""]
    for i, s in enumerate(v.columns):
        mk = markers[blevels.index(batches.get(s, "")) % len(markers)] if batches else "o"
        ax.scatter(pcs[i, 0], pcs[i, 1], s=70, color=colors[groups[s]], marker=mk, edgecolor="white", linewidth=0.8,
                   zorder=3)
        ax.annotate(s, (pcs[i, 0], pcs[i, 1]), fontsize=7, xytext=(4, 4), textcoords="offset points", color="#3a4350")
    handles = [matplotlib.lines.Line2D([], [], marker="o", ls="", color=c, label=g) for g, c in colors.items()]
    if batches and len(blevels) > 1:
        handles += [matplotlib.lines.Line2D([], [], marker=markers[i % len(markers)], ls="", color="#666",
                                            label=f"batch {b}") for i, b in enumerate(blevels)]
    ax.legend(handles=handles, frameon=False, fontsize=8, loc="best")
    ax.set_xlabel(f"PC1 ({var[0] * 100:.1f}% variance)")
    ax.set_ylabel(f"PC2 ({var[1] * 100:.1f}% variance)" if len(var) > 1 else "PC2")
    ax.set_title(f"PCA on the {len(v)} most variable genes (VST)", fontsize=11)
    ax.axhline(0, color="#e3e6ea", lw=0.8, zorder=0)
    ax.axvline(0, color="#e3e6ea", lw=0.8, zorder=0)
    _style(ax)


def plot_volcano(fig, res, alpha, thr, lfc, title, n_label=15):
    ax = fig.add_subplot(111)
    d = res.dropna(subset=["padj", lfc])
    x, y = d[lfc], -np.log10(d["padj"].clip(lower=1e-300))
    up = (d["padj"] < alpha) & (x >= thr)
    down = (d["padj"] < alpha) & (x <= -thr)
    ns = ~(up | down)
    ax.scatter(x[ns], y[ns], s=5, color=NS_COLOR, alpha=0.6, linewidths=0, rasterized=True)
    ax.scatter(x[up], y[up], s=9, color=UP_COLOR, alpha=0.85, linewidths=0, label=f"UP ({int(up.sum())})",
               rasterized=True)
    ax.scatter(x[down], y[down], s=9, color=DOWN_COLOR, alpha=0.85, linewidths=0, label=f"DOWN ({int(down.sum())})",
               rasterized=True)
    ax.axhline(-np.log10(alpha), color="#8a929e", lw=0.8, ls="--")
    if thr > 0:
        for t in (-thr, thr):
            ax.axvline(t, color="#8a929e", lw=0.8, ls="--")
    top = d[up | down].nsmallest(n_label, "padj")
    for gid, row in top.iterrows():
        ax.annotate(str(row.get("symbol", gid)), (row[lfc], -np.log10(max(row["padj"], 1e-300))), fontsize=7,
                    xytext=(3, 3), textcoords="offset points", color="#222")
    ax.set_xlabel("log2 fold change" + (" (shrunk)" if "shrunk" in lfc else ""))
    ax.set_ylabel("−log10 padj")
    ax.set_title(title, fontsize=11)
    ax.legend(frameon=False, fontsize=8, loc="upper left")
    _style(ax)


def plot_ma(fig, res, alpha, thr, lfc):
    ax = fig.add_subplot(111)
    d = res.dropna(subset=[lfc])
    d = d[d["baseMean"] > 0]
    sig = (d["padj"] < alpha) & (d[lfc].abs() >= thr)
    ax.scatter(d.loc[~sig, "baseMean"], d.loc[~sig, lfc], s=5, color=NS_COLOR, alpha=0.6, linewidths=0, rasterized=True)
    ax.scatter(d.loc[sig & (d[lfc] > 0), "baseMean"], d.loc[sig & (d[lfc] > 0), lfc], s=8, color=UP_COLOR,
               linewidths=0, rasterized=True)
    ax.scatter(d.loc[sig & (d[lfc] < 0), "baseMean"], d.loc[sig & (d[lfc] < 0), lfc], s=8, color=DOWN_COLOR,
               linewidths=0, rasterized=True)
    ax.set_xscale("log")
    ax.axhline(0, color="#444", lw=0.8)
    ax.set_xlabel("Mean normalised expression (baseMean)")
    ax.set_ylabel("log2 fold change")
    ax.set_title("MA plot", fontsize=11)
    _style(ax)


def plot_heatmap(fig, vst, res, groups, ref, case, n, alpha):
    from scipy.cluster.hierarchy import linkage, leaves_list
    sig = res[(res["padj"] < alpha) & (res["direction"] != "")]
    pick = (sig if len(sig) >= 2 else res.dropna(subset=["padj"])).nsmallest(n, "padj")
    genes = [g for g in pick.index if g in vst.index]
    cols = [s for s in _ordered_samples(groups, ref, case) if s in vst.columns]
    m = vst.loc[genes, cols]
    z = m.sub(m.mean(axis=1), axis=0).div(m.std(axis=1).replace(0, 1), axis=0)
    if len(z) > 2:
        z = z.iloc[leaves_list(linkage(z.values, "average", metric="euclidean"))]
    gs = fig.add_gridspec(2, 2, height_ratios=[1, 28], width_ratios=[30, 1])
    axg, ax, cax = fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[1, 0]), fig.add_subplot(gs[1, 1])
    colors = _group_colors(groups, ref, case)
    glist = list(colors)
    axg.imshow([[glist.index(groups[s]) for s in cols]], aspect="auto",
               cmap=ListedColormap([colors[g] for g in glist]), vmin=-0.5, vmax=len(glist) - 0.5)
    axg.set_xticks([])
    axg.set_yticks([])
    from matplotlib.patches import Patch
    axg.legend(handles=[Patch(color=colors[g], label=g) for g in glist], ncol=min(len(glist), 6), loc="lower left",
               bbox_to_anchor=(0, 1.0), frameon=False, fontsize=8, handlelength=1, borderaxespad=0.2)
    im = ax.imshow(z.values, aspect="auto", cmap="RdBu_r", vmin=-2.5, vmax=2.5, interpolation="nearest")
    ax.set_xticks(range(len(cols)))
    ax.set_xticklabels(cols, rotation=60, ha="right", fontsize=7)
    labels = res["symbol"].reindex(z.index).fillna(pd.Series(z.index, index=z.index)) if "symbol" in res else z.index
    ax.set_yticks(range(len(z)))
    ax.set_yticklabels(labels, fontsize=6 if len(z) > 40 else 7)
    fig.colorbar(im, cax=cax, label="z-score (VST)")


def plot_corr(fig, vst, groups, ref, case):
    ax = fig.add_subplot(111)
    cols = [s for s in _ordered_samples(groups, ref, case) if s in vst.columns]
    c = vst[cols].corr(method="pearson")
    im = ax.imshow(c.values, cmap="viridis", vmin=max(0.0, float(np.nanmin(c.values))), vmax=1)
    ax.set_xticks(range(len(cols)))
    ax.set_xticklabels(cols, rotation=60, ha="right", fontsize=7)
    ax.set_yticks(range(len(cols)))
    ax.set_yticklabels(cols, fontsize=7)
    if len(cols) <= 16:
        for i in range(len(cols)):
            for j in range(len(cols)):
                ax.text(j, i, f"{c.values[i, j]:.2f}", ha="center", va="center", fontsize=6, color="white")
    fig.colorbar(im, ax=ax, label="Pearson r (VST)")
    ax.set_title("Sample correlation", fontsize=11)


def _short(t, n=60):
    t = str(t)
    return t if len(t) <= n else t[: n - 1] + "…"


def plot_ora_bars(fig, ora_df, n=10):
    dirs = [d for d in ("UP", "DOWN") if (ora_df["Direction"] == d).any()] if len(ora_df) else []
    if not dirs:
        fig.text(0.5, 0.5, "No ORA results", ha="center", va="center")
        return
    for i, d in enumerate(dirs):
        ax = fig.add_subplot(len(dirs), 1, i + 1)
        sub = ora_df[ora_df["Direction"] == d].nsmallest(n, "FDR").iloc[::-1]
        y = -np.log10(sub["FDR"].clip(lower=1e-300))
        ax.barh(range(len(sub)), y, color=UP_COLOR if d == "UP" else DOWN_COLOR, height=0.7)
        ax.set_yticks(range(len(sub)))
        ax.set_yticklabels([f"{_short(t)}  [{lib}]" for t, lib in zip(sub["Term"], sub["Library"])], fontsize=7)
        ax.axvline(-np.log10(0.05), color="#8a929e", lw=0.8, ls="--")
        ax.set_xlabel("−log10 FDR")
        ax.set_title(f"Enriched terms in {d}", fontsize=10, loc="left")
        _style(ax)


def plot_gsea_bars(fig, gsea_df, n=10, fdr=0.25):
    ax = fig.add_subplot(111)
    d = gsea_df.dropna(subset=["NES"])
    pos = d[d["NES"] > 0].nsmallest(n, "FDR")
    neg = d[d["NES"] < 0].nsmallest(n, "FDR")
    sub = pd.concat([neg, pos]).sort_values("NES")
    if sub.empty:
        fig.text(0.5, 0.5, "No GSEA results", ha="center", va="center")
        return
    colors = [UP_COLOR if v > 0 else DOWN_COLOR for v in sub["NES"]]
    alphas = [1.0 if f < fdr else 0.35 for f in sub["FDR"]]
    for i, (v, c, a) in enumerate(zip(sub["NES"], colors, alphas)):
        ax.barh(i, v, color=c, alpha=a, height=0.7)
    ax.set_yticks(range(len(sub)))
    ax.set_yticklabels([f"{_short(t)}  [{lib}]" for t, lib in zip(sub["Term"], sub["Library"])], fontsize=7)
    ax.axvline(0, color="#444", lw=0.8)
    ax.set_xlabel("NES (solid colour: FDR < %.2f)" % fdr)
    ax.set_title("Preranked GSEA: top terms by direction", fontsize=10, loc="left")
    _style(ax)


def plot_enrichment(fig, ranking, info, term, lib):
    hits = info.get("hits") or info.get("hit_indices")
    res_curve = info.get("RES")
    if hits is None or res_curve is None:
        raise ValueError("The installed gseapy version does not expose the running ES for this term.")
    res_curve = np.asarray(res_curve, float)
    gs = fig.add_gridspec(3, 1, height_ratios=[3, 0.6, 1.6], hspace=0.05)
    a1, a2, a3 = fig.add_subplot(gs[0]), fig.add_subplot(gs[1]), fig.add_subplot(gs[2])
    x = np.arange(len(res_curve))
    a1.plot(x, res_curve, color="#1f8a70", lw=1.6)
    a1.axhline(0, color="#999", lw=0.7)
    a1.set_ylabel("Running ES")
    a1.set_xticks([])
    nes, fdr = info.get("nes", np.nan), info.get("fdr", np.nan)
    a1.set_title(f"{_short(term, 80)}\n{lib} — NES {float(nes):.2f}, FDR {float(fdr):.3g}", fontsize=10, loc="left")
    a2.vlines(hits, 0, 1, color="#222", lw=0.6)
    a2.set_yticks([])
    a2.set_xticks([])
    a2.set_xlim(a1.get_xlim())
    a3.fill_between(np.arange(len(ranking)), ranking.values, color="#9aa3af", lw=0)
    a3.axhline(0, color="#444", lw=0.6)
    a3.set_xlim(a1.get_xlim())
    a3.set_ylabel("Metric")
    a3.set_xlabel("Position in the ranked gene list")
    for a in (a1, a3):
        _style(a)


# =============================================================================
# 5. GUI components
# =============================================================================
QSS = """
* { font-family: 'IBM Plex Sans', 'Source Sans 3', 'Segoe UI', 'Helvetica Neue', Arial, sans-serif; font-size: 13px; }
QMainWindow, QWidget { background: #eef1f4; color: #1e2a36; }
#header { background: #ffffff; border-bottom: 1px solid #d9dee5; }
#appTitle { font-size: 19px; font-weight: 700; color: #1e2a36; }
#appSub { color: #6b7684; }
#nav { background: #1e2a36; border: none; padding: 14px 8px; outline: 0; }
#nav::item { color: #b7c2cf; padding: 12px 14px; border-radius: 6px; margin: 2px 2px; }
#nav::item:selected { background: #4b3fa8; color: #ffffff; }
#nav::item:hover:!selected { background: #2a3847; }
#pageTitle { font-size: 21px; font-weight: 600; color: #1e2a36; }
#pageSub { color: #5d6875; }
QGroupBox { background: #ffffff; border: 1px solid #d9dee5; border-radius: 8px; margin-top: 16px; padding: 14px 12px 10px 12px; }
QGroupBox::title { subcontrol-origin: margin; left: 10px; padding: 0 4px; color: #4b3fa8; font-weight: 600; }
QLineEdit, QSpinBox, QDoubleSpinBox, QComboBox { background: #ffffff; border: 1px solid #c9d0d9; border-radius: 5px; padding: 4px 6px; min-height: 22px; }
QLineEdit:focus, QSpinBox:focus, QDoubleSpinBox:focus, QComboBox:focus { border: 1px solid #4b3fa8; }
QPushButton, QToolButton { background: #ffffff; border: 1px solid #c9d0d9; border-radius: 5px; padding: 6px 12px; color: #1e2a36; }
QPushButton:hover, QToolButton:hover { background: #f2f0fb; border-color: #4b3fa8; }
QPushButton:disabled { color: #a0a8b3; background: #f4f6f8; }
QPushButton#primary { background: #4b3fa8; border: none; color: #ffffff; font-weight: 600; padding: 9px 20px; }
QPushButton#primary:hover { background: #5b4fc0; }
QPushButton#primary:disabled { background: #a9a3d6; }
QTableWidget, QTableView, QListWidget { background: #ffffff; alternate-background-color: #f6f7f9; border: 1px solid #d9dee5; border-radius: 6px; gridline-color: #e7eaee; selection-background-color: #ded9f6; selection-color: #1e2a36; }
QHeaderView::section { background: #f2f4f7; color: #3b4652; border: none; border-right: 1px solid #e1e5ea; border-bottom: 1px solid #d9dee5; padding: 5px 6px; font-weight: 600; }
QTabWidget::pane { border: 1px solid #d9dee5; border-radius: 6px; background: #ffffff; top: -1px; }
QTabBar::tab { background: transparent; color: #6b7684; padding: 8px 16px; border-bottom: 2px solid transparent; }
QTabBar::tab:selected { color: #4b3fa8; border-bottom: 2px solid #4b3fa8; }
QProgressBar { background: #ffffff; border: 1px solid #c9d0d9; border-radius: 5px; text-align: center; min-height: 16px; }
QProgressBar::chunk { background: #4b3fa8; border-radius: 4px; }
#logPanel { background: #ffffff; border-top: 1px solid #d9dee5; }
QPlainTextEdit#log { background: #f8f9fb; border: 1px solid #e1e5ea; border-radius: 6px; font-family: 'JetBrains Mono', 'Cascadia Mono', 'DejaVu Sans Mono', 'Menlo', monospace; font-size: 12px; color: #2c3640; }
QScrollArea { border: none; background: transparent; }
QLabel, QCheckBox, QRadioButton { background: transparent; }
QCheckBox::indicator, QRadioButton::indicator { width: 14px; height: 14px; border: 1px solid #9aa3af; background: #ffffff; }
QCheckBox::indicator { border-radius: 3px; }
QRadioButton::indicator { border-radius: 8px; }
QCheckBox::indicator:checked, QRadioButton::indicator:checked { background: #4b3fa8; border: 1px solid #4b3fa8; }
QSplitter::handle { background: #d9dee5; }
QToolTip { background: #1e2a36; color: #ffffff; border: none; padding: 6px; }
#hint { color: #6b7684; font-size: 12px; }
QPushButton#role_ctrl { color: #1f5a99; background: #e4eef9; border-color: #a9c4e8; font-weight: 600; }
QPushButton#role_case { color: #b02a3f; background: #fbe7ea; border-color: #ebb2bc; font-weight: 600; }
QPushButton#role_other { color: #4f5b68; }
QPushButton#role_excl { color: #8a929e; background: #f3f4f6; }
QToolButton#chip_ctrl, QToolButton#chip_case, QToolButton#chip_other, QToolButton#chip_excl, QToolButton#chip_mixed { border-radius: 12px; padding: 4px 12px; }
QToolButton#chip_ctrl { background: #e4eef9; color: #1f5a99; border-color: #a9c4e8; }
QToolButton#chip_case { background: #fbe7ea; color: #b02a3f; border-color: #ebb2bc; }
QToolButton#chip_excl { background: #f3f4f6; color: #8a929e; }
#designSummary { padding: 6px 8px; background: #f6f7f9; border: 1px solid #e1e5ea; border-radius: 6px; }
#summary { font-size: 14px; font-weight: 600; color: #1e2a36; padding: 4px 0; }
"""


def apply_palette(app):
    pal = QPalette()
    pal.setColor(QPalette.ColorRole.Window, QColor("#eef1f4"))
    pal.setColor(QPalette.ColorRole.Base, QColor("#ffffff"))
    pal.setColor(QPalette.ColorRole.AlternateBase, QColor("#f6f7f9"))
    pal.setColor(QPalette.ColorRole.Text, QColor("#1e2a36"))
    pal.setColor(QPalette.ColorRole.WindowText, QColor("#1e2a36"))
    pal.setColor(QPalette.ColorRole.Button, QColor("#ffffff"))
    pal.setColor(QPalette.ColorRole.ButtonText, QColor("#1e2a36"))
    pal.setColor(QPalette.ColorRole.Highlight, QColor("#4b3fa8"))
    pal.setColor(QPalette.ColorRole.HighlightedText, QColor("#ffffff"))
    app.setPalette(pal)


def fmt_value(v):
    if isinstance(v, (float, np.floating)):
        if np.isnan(v):
            return ""
        if v != 0 and abs(v) < 1e-3:
            return f"{v:.2e}"
        return f"{v:.4g}" if abs(v) < 1e5 else f"{v:,.0f}"
    if isinstance(v, (bool, np.bool_)):
        return "yes" if v else ""
    return str(v)


class PandasModel(QAbstractTableModel):
    def __init__(self, df=None):
        super().__init__()
        self._df = df if df is not None else pd.DataFrame()

    def set_df(self, df):
        self.beginResetModel()
        self._df = df.reset_index(drop=True)
        self.endResetModel()

    def rowCount(self, parent=QModelIndex()):
        return 0 if parent.isValid() else len(self._df)

    def columnCount(self, parent=QModelIndex()):
        return 0 if parent.isValid() else self._df.shape[1]

    def data(self, index, role=Qt.ItemDataRole.DisplayRole):
        if not index.isValid():
            return None
        v = self._df.iat[index.row(), index.column()]
        if role == Qt.ItemDataRole.DisplayRole:
            return fmt_value(v)
        if role == Qt.ItemDataRole.UserRole:
            if isinstance(v, (bool, np.bool_)):
                return int(v)
            if isinstance(v, (int, np.integer)):
                return int(v)
            if isinstance(v, (float, np.floating)):
                return float("inf") if np.isnan(v) else float(v)
            return str(v)
        if role == Qt.ItemDataRole.TextAlignmentRole and isinstance(v, (int, float, np.number)):
            return int(Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
        if role == Qt.ItemDataRole.ForegroundRole and v in ("UP", "DOWN"):
            return QBrush(QColor(UP_COLOR if v == "UP" else DOWN_COLOR))
        return None

    def headerData(self, section, orientation, role=Qt.ItemDataRole.DisplayRole):
        if role != Qt.ItemDataRole.DisplayRole:
            return None
        if orientation == Qt.Orientation.Horizontal:
            return str(self._df.columns[section])
        return str(section + 1)


def make_table():
    view, model, proxy = QTableView(), PandasModel(), QSortFilterProxyModel()
    proxy.setSourceModel(model)
    proxy.setSortRole(Qt.ItemDataRole.UserRole)
    proxy.setFilterCaseSensitivity(Qt.CaseSensitivity.CaseInsensitive)
    proxy.setFilterKeyColumn(-1)
    view.setModel(proxy)
    view.setSortingEnabled(True)
    view.setAlternatingRowColors(True)
    view.setSelectionBehavior(QAbstractItemView.SelectionBehavior.SelectRows)
    view.horizontalHeader().setStretchLastSection(True)
    view.verticalHeader().setDefaultSectionSize(24)
    view.verticalHeader().setVisible(False)
    return view, model, proxy


class PathEdit(QWidget):
    def __init__(self, mode="file", filt="All files (*)", placeholder=""):
        super().__init__()
        lay = QHBoxLayout(self)
        lay.setContentsMargins(0, 0, 0, 0)
        lay.setSpacing(4)
        self.le = QLineEdit()
        self.le.setPlaceholderText(placeholder)
        btn = QToolButton()
        btn.setText("Browse…")
        btn.clicked.connect(self.browse)
        lay.addWidget(self.le, 1)
        lay.addWidget(btn)
        self.mode, self.filt = mode, filt

    def browse(self):
        start = self.le.text() or str(Path.home())
        if self.mode == "dir":
            p = QFileDialog.getExistingDirectory(self, "Select folder", start)
        else:
            p, _ = QFileDialog.getOpenFileName(self, "Select file", str(Path(start).parent) if
                                               Path(start).is_file() else start, self.filt)
        if p:
            self.le.setText(p)

    def text(self):
        return self.le.text().strip()

    def setText(self, t):
        self.le.setText(str(t))


class PlotPanel(QWidget):
    def __init__(self):
        super().__init__()
        lay = QVBoxLayout(self)
        top = QHBoxLayout()
        self.combo = QComboBox()
        self.combo.currentTextChanged.connect(self.draw)
        self.save_btn = QPushButton("Export all (PNG 300 dpi + PDF)")
        self.save_btn.clicked.connect(self._export_clicked)
        top.addWidget(QLabel("Plot"))
        top.addWidget(self.combo, 1)
        top.addWidget(self.save_btn)
        self.fig = Figure(figsize=(8, 6), dpi=100, layout="constrained", facecolor="white")
        self.canvas = FigureCanvas(self.fig)
        self.toolbar = NavToolbar(self.canvas, self)
        lay.addLayout(top)
        lay.addWidget(self.toolbar)
        lay.addWidget(self.canvas, 1)
        self.plots, self.export_dir = {}, None

    def set_plots(self, plots, export_dir=None):
        self.plots, self.export_dir = plots, export_dir
        self.combo.blockSignals(True)
        self.combo.clear()
        self.combo.addItems(list(plots))
        self.combo.blockSignals(False)
        if plots:
            self.draw(self.combo.currentText())

    def draw(self, name):
        if name not in self.plots:
            return
        self.fig.clear()
        try:
            self.plots[name](self.fig)
        except Exception as e:
            self.fig.clear()
            self.fig.text(0.5, 0.5, f"Plot not available:\n{e}", ha="center", va="center", wrap=True)
        self.canvas.draw_idle()

    def export_all(self, outdir):
        return export_plots(self.plots, outdir)

    def _export_clicked(self):
        d = QFileDialog.getExistingDirectory(self, "Folder for the plots", self.export_dir or str(Path.home()))
        if d:
            n = self.export_all(d)
            QMessageBox.information(self, APP_NAME, f"{n} plots saved in {d}")


def card(title):
    box = QGroupBox(title)
    return box


def hint(text):
    lab = QLabel(text)
    lab.setObjectName("hint")
    lab.setWordWrap(True)
    return lab


def scroll(widget):
    sa = QScrollArea()
    sa.setWidgetResizable(True)
    sa.setWidget(widget)
    sa.setFrameShape(QFrame.Shape.NoFrame)
    sa.setHorizontalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAlwaysOff)
    widget.setMinimumWidth(widget.sizeHint().width())
    return sa


class Page(QWidget):
    def __init__(self, mw, title, subtitle):
        super().__init__()
        self.mw = mw
        lay = QVBoxLayout(self)
        lay.setContentsMargins(20, 14, 20, 8)
        t = QLabel(title)
        t.setObjectName("pageTitle")
        s = QLabel(subtitle)
        s.setObjectName("pageSub")
        s.setWordWrap(True)
        lay.addWidget(t)
        lay.addWidget(s)
        self.body = QVBoxLayout()
        lay.addLayout(self.body, 1)


# ----------------------------------------------------------------------------- Environment
class EnvPage(Page):
    def __init__(self, mw):
        super().__init__(mw, "Environment", "Dependency status. The Differential expression and Pathway steps "
                                         "work with the Python libraries alone; alignment and counting require "
                                         "the external tools, which can be installed here with one click (Linux/macOS).")
        row = QHBoxLayout()
        b1 = card("Python libraries")
        l1 = QVBoxLayout(b1)
        self.py_table = QTableWidget(0, 2)
        self.py_table.setHorizontalHeaderLabels(["Package", "Version"])
        self.py_table.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeMode.Stretch)
        self.py_table.verticalHeader().setVisible(False)
        l1.addWidget(self.py_table)
        b2 = card("External tools")
        l2 = QVBoxLayout(b2)
        self.tool_table = QTableWidget(0, 2)
        self.tool_table.setHorizontalHeaderLabels(["Tool", "Path"])
        self.tool_table.horizontalHeader().setSectionResizeMode(1, QHeaderView.ResizeMode.Stretch)
        self.tool_table.verticalHeader().setVisible(False)
        l2.addWidget(self.tool_table)
        btns = QHBoxLayout()
        self.install_btn = QPushButton("Install / update external tools")
        self.install_btn.setObjectName("primary")
        self.install_btn.clicked.connect(self.install)
        rb = QPushButton("Check again")
        rb.clicked.connect(self.refresh)
        btns.addWidget(self.install_btn)
        btns.addWidget(rb)
        btns.addStretch()
        l2.addLayout(btns)
        l2.addWidget(hint(f"Tools are installed from Bioconda into {TOOLS_ENV}, without administrator rights, "
                          f"with commands available in {LOCAL_BIN}. Tools already on the PATH are used as they are."))
        row.addWidget(b1, 2)
        row.addWidget(b2, 3)
        self.body.addLayout(row)
        self.body.addWidget(hint(f"Application folder: {APP_HOME} — Python environment: {sys.prefix} — "
                                 f"Python {platform.python_version()} on {platform.system()} {platform.machine()}"))
        self.body.addStretch()
        self.refresh()

    def refresh(self):
        pk = {"PySide6": "PySide6", "pandas": "pandas", "numpy": "numpy", "scipy": "scipy",
              "matplotlib": "matplotlib", "pydeseq2": "pydeseq2", "gseapy": "gseapy"}
        self.py_table.setRowCount(len(pk))
        for i, (m, dist) in enumerate(pk.items()):
            self.py_table.setItem(i, 0, QTableWidgetItem(m))
            self.py_table.setItem(i, 1, QTableWidgetItem(_pkg_version(dist) or "not installed"))
        self.tool_table.setRowCount(len(EXT_TOOLS))
        for i, t in enumerate(EXT_TOOLS):
            p = tool_path(t)
            self.tool_table.setItem(i, 0, QTableWidgetItem(t))
            it = QTableWidgetItem(p or "not found")
            it.setForeground(QBrush(QColor("#1f8a70" if p else UP_COLOR)))
            self.tool_table.setItem(i, 1, it)

    def install(self):
        self.mw.start_task(task_install_tools, title="External tool installation",
                           on_done=lambda r: (self.refresh(),
                                              QMessageBox.information(self, APP_NAME, "External tools installed.")))


# ----------------------------------------------------------------------------- Alignment
class AlignPage(Page):
    def __init__(self, mw):
        super().__init__(mw, "1  Alignment", "Detects FASTQ files (multi-lane, paired- or single-end), "
                                               "runs optional trimming and aligns with HISAT2 or STAR. "
                                               "BAM files already present are skipped, so an interrupted analysis "
                                               "resumes where it stopped.")
        split = QSplitter(Qt.Orientation.Horizontal)
        left = card("FASTQ samples")
        ll = QVBoxLayout(left)
        r = QHBoxLayout()
        self.fq_dir = PathEdit("dir", placeholder="Folder containing the FASTQ files")
        det = QPushButton("Detect samples")
        det.clicked.connect(self.detect)
        r.addWidget(self.fq_dir, 1)
        r.addWidget(det)
        ll.addLayout(r)
        self.table = QTableWidget(0, 5)
        self.table.setHorizontalHeaderLabels(["Use", "Sample", "Type", "R1", "R2"])
        hh = self.table.horizontalHeader()
        hh.setSectionResizeMode(3, QHeaderView.ResizeMode.Stretch)
        hh.setSectionResizeMode(4, QHeaderView.ResizeMode.Stretch)
        self.table.verticalHeader().setVisible(False)
        ll.addWidget(self.table, 1)
        ll.addWidget(hint("The sample name is editable and becomes the name of the BAM file and of the column "
                          "in the count matrix."))
        split.addWidget(left)

        right = QWidget()
        rl = QVBoxLayout(right)
        rl.setContentsMargins(0, 0, 0, 0)
        ref = card("Reference")
        f = QFormLayout(ref)
        self.aligner = QComboBox()
        self.aligner.addItems(["HISAT2", "STAR"])
        self.aligner.setToolTip("HISAT2: little RAM (~8 GB). STAR: faster and more sensitive on junctions, "
                                "~32 GB of RAM for mammalian genomes.")
        self.fasta = PathEdit(filt="FASTA (*.fa *.fasta *.fna *.fa.gz *.fasta.gz *.fna.gz);;All (*)",
                              placeholder="Genome .fa / .fa.gz")
        self.gtf = PathEdit(filt="GTF (*.gtf *.gtf.gz);;All (*)", placeholder="Annotation .gtf / .gtf.gz")
        self.index = PathEdit(placeholder="Optional: existing index (.1.ht2 file or STAR folder)")
        self.index_dir_btn = QToolButton()
        self.index_dir_btn.setText("Folder…")
        self.index_dir_btn.clicked.connect(self._pick_index_dir)
        irow = QHBoxLayout()
        irow.addWidget(self.index, 1)
        irow.addWidget(self.index_dir_btn)
        self.overhang = QSpinBox()
        self.overhang.setRange(25, 500)
        self.overhang.setValue(100)
        self.overhang.setToolTip("STAR --sjdbOverhang: read length − 1 (100 works in most cases).")
        f.addRow("Aligner", self.aligner)
        f.addRow("Genome FASTA", self.fasta)
        f.addRow("GTF annotation", self.gtf)
        f.addRow("Index", irow)
        f.addRow("sjdbOverhang (STAR)", self.overhang)
        rl.addWidget(ref)
        opt = card("Options")
        fo = QFormLayout(opt)
        self.trim = QCheckBox("Adapter and quality trimming with fastp")
        self.trim.setChecked(True)
        self.extra = QLineEdit()
        self.extra.setPlaceholderText("Extra arguments for the aligner, e.g. --rna-strandness RF")
        self.skip_err = QCheckBox("If a sample fails, skip it and continue with the others")
        self.skip_err.setChecked(True)
        fo.addRow(self.trim)
        fo.addRow(self.skip_err)
        fo.addRow("Extra arguments", self.extra)
        fo.addRow(hint("Completed steps (trimming, index, BAM) are recognised and skipped: "
                       "when relaunched, the analysis restarts from the first unfinished sample."))
        rl.addWidget(opt)
        run = QPushButton("Align selected samples")
        run.setObjectName("primary")
        run.clicked.connect(self.run)
        rl.addWidget(run)
        rl.addStretch()
        split.addWidget(scroll(right))
        split.setSizes([700, 480])
        self.body.addWidget(split, 1)

    def _pick_index_dir(self):
        d = QFileDialog.getExistingDirectory(self, "STAR index folder", self.index.text() or str(Path.home()))
        if d:
            self.index.setText(d)

    def detect(self):
        d = self.fq_dir.text()
        if not d or not Path(d).is_dir():
            QMessageBox.warning(self, APP_NAME, "Select a valid folder.")
            return
        samples = discover_fastq(d)
        self.table.setRowCount(0)
        for name, fq in samples.items():
            r = self.table.rowCount()
            self.table.insertRow(r)
            use = QTableWidgetItem()
            use.setFlags(Qt.ItemFlag.ItemIsUserCheckable | Qt.ItemFlag.ItemIsEnabled)
            use.setCheckState(Qt.CheckState.Checked)
            self.table.setItem(r, 0, use)
            it = QTableWidgetItem(name)
            it.setData(Qt.ItemDataRole.UserRole, fq)
            self.table.setItem(r, 1, it)
            for c, txt in ((2, "PE" if fq["R2"] else "SE"),
                           (3, ", ".join(Path(p).name for p in fq["R1"])),
                           (4, ", ".join(Path(p).name for p in fq["R2"]))):
                x = QTableWidgetItem(txt)
                x.setFlags(Qt.ItemFlag.ItemIsEnabled | Qt.ItemFlag.ItemIsSelectable)
                self.table.setItem(r, c, x)
        self.table.resizeColumnToContents(1)
        self.mw.log(f"Detected {len(samples)} samples in {d}")
        if not samples:
            QMessageBox.information(self, APP_NAME, "No FASTQ found (*.fastq, *.fq, also .gz).")

    def run(self):
        samples = {}
        for r in range(self.table.rowCount()):
            if self.table.item(r, 0).checkState() != Qt.CheckState.Checked:
                continue
            name = sanitize(self.table.item(r, 1).text())
            if name in samples:
                QMessageBox.warning(self, APP_NAME, f"Duplicate sample name: {name}")
                return
            samples[name] = self.table.item(r, 1).data(Qt.ItemDataRole.UserRole)
        if not samples:
            QMessageBox.warning(self, APP_NAME, "No sample selected.")
            return
        cfg = dict(outdir=self.mw.outdir(), threads=self.mw.threads(), aligner=self.aligner.currentText(),
                   fasta=self.fasta.text(), gtf=self.gtf.text(), index=self.index.text(),
                   overhang=self.overhang.value(), trim=self.trim.isChecked(), extra=self.extra.text(),
                   skip_errors=self.skip_err.isChecked())
        self.mw.start_task(task_align, samples, cfg, title="Alignment", on_done=self.done)

    def done(self, r):
        self.mw.count_page.receive_bams(r["bams"], r["paired"], r["gtf"])
        self.mw.goto(self.mw.count_page)
        msg = f"Alignment completed: {len(r['bams'])} BAM files, loaded into the Counts step."
        if r.get("failed"):
            msg += "\n\nSamples not aligned (excluded):\n" + "\n".join(f"• {k}: {v[:120]}"
                                                                         for k, v in r["failed"].items())
            QMessageBox.warning(self, APP_NAME, msg)
        else:
            QMessageBox.information(self, APP_NAME, msg)


# ----------------------------------------------------------------------------- Counts
class CountPage(Page):
    def __init__(self, mw):
        super().__init__(mw, "2  Counts", "Counts reads per gene with featureCounts starting from BAM files "
                                           "(produced here or elsewhere). Produces a raw count matrix, TPM, "
                                           "gene annotation and an assignment summary.")
        split = QSplitter(Qt.Orientation.Horizontal)
        left = card("BAM files")
        ll = QVBoxLayout(left)
        self.bam_list = QListWidget()
        self.bam_list.setSelectionMode(QAbstractItemView.SelectionMode.ExtendedSelection)
        ll.addWidget(self.bam_list, 1)
        br = QHBoxLayout()
        for txt, fn in (("Add BAM…", self.add_files), ("Add folder…", self.add_dir),
                        ("Remove", self.remove), ("Clear", self.bam_list.clear)):
            b = QPushButton(txt)
            b.clicked.connect(fn)
            br.addWidget(b)
        ll.addLayout(br)
        split.addWidget(left)

        right = QWidget()
        rl = QVBoxLayout(right)
        rl.setContentsMargins(0, 0, 0, 0)
        par = card("featureCounts parameters")
        f = QFormLayout(par)
        self.gtf = PathEdit(filt="GTF (*.gtf *.gtf.gz);;All (*)", placeholder="Annotation .gtf")
        self.strand = QComboBox()
        self.strand.addItems(["0 · non stranded", "1 · stranded (forward)",
                              "2 · reverse stranded (TruSeq/Illumina Stranded mRNA, NEBNext Directional)"])
        self.strand.setCurrentIndex(2)
        self.strand.setToolTip("If the % of assigned reads is very low (<40%), try another option.")
        self.paired = QComboBox()
        self.paired.addItems(["Detect automatically", "Paired-end", "Single-end"])
        self.feature = QLineEdit("exon")
        self.attribute = QLineEdit("gene_id")
        self.mapq = QSpinBox()
        self.mapq.setRange(0, 255)
        self.mapq.setValue(10)
        self.mapq.setToolTip("Minimum mapping quality. A value of 10 effectively excludes multi-mapping reads "
                             "(HISAT2 = 1, STAR = 0/1/3).")
        f.addRow("GTF", self.gtf)
        f.addRow("Strandedness (-s)", self.strand)
        f.addRow("Library", self.paired)
        f.addRow("Feature (-t)", self.feature)
        f.addRow("Attribute (-g)", self.attribute)
        f.addRow("Minimum MAPQ (-Q)", self.mapq)
        self.skip_bad = QCheckBox("Automatically exclude damaged BAM files and continue")
        self.skip_bad.setChecked(True)
        f.addRow(self.skip_bad)
        rl.addWidget(par)
        run = QPushButton("Compute counts")
        run.setObjectName("primary")
        run.clicked.connect(self.run)
        rl.addWidget(run)
        res = card("Assignment summary")
        rsl = QVBoxLayout(res)
        self.sum_view, self.sum_model, _ = make_table()
        rsl.addWidget(self.sum_view)
        self.to_de = QPushButton("Use these counts for differential expression")
        self.to_de.setEnabled(False)
        self.to_de.clicked.connect(lambda: self.mw.goto(self.mw.de_page))
        rsl.addWidget(self.to_de)
        rl.addWidget(res, 1)
        split.addWidget(scroll(right))
        split.setSizes([520, 660])
        self.body.addWidget(split, 1)

    def receive_bams(self, bams, paired, gtf):
        self.bam_list.clear()
        self.bam_list.addItems(bams)
        self.paired.setCurrentIndex(1 if paired else 2)
        if gtf:
            self.gtf.setText(gtf)

    def add_files(self):
        fs, _ = QFileDialog.getOpenFileNames(self, "Select BAM files", str(Path.home()), "BAM (*.bam)")
        self.bam_list.addItems([f for f in fs if not self.bam_list.findItems(f, Qt.MatchFlag.MatchExactly)])

    def add_dir(self):
        d = QFileDialog.getExistingDirectory(self, "Folder containing BAM files", str(Path.home()))
        if d:
            fs = sorted(str(p) for p in Path(d).rglob("*.bam"))
            self.bam_list.addItems([f for f in fs if not self.bam_list.findItems(f, Qt.MatchFlag.MatchExactly)])

    def remove(self):
        for it in self.bam_list.selectedItems():
            self.bam_list.takeItem(self.bam_list.row(it))

    def run(self):
        bams = [self.bam_list.item(i).text() for i in range(self.bam_list.count())]
        if not bams:
            QMessageBox.warning(self, APP_NAME, "Add at least one BAM file.")
            return
        cfg = dict(outdir=self.mw.outdir(), threads=self.mw.threads(), gtf=self.gtf.text(),
                   strand=self.strand.currentIndex(), paired={0: None, 1: True, 2: False}[self.paired.currentIndex()],
                   feature=self.feature.text() or "exon", attribute=self.attribute.text() or "gene_id",
                   mapq=self.mapq.value(), skip_bad=self.skip_bad.isChecked())
        self.mw.start_task(task_featurecounts, bams, cfg, title="featureCounts counting", on_done=self.done)

    def done(self, r):
        if r["summary"] is not None:
            self.sum_model.set_df(r["summary"])
            self.sum_view.resizeColumnsToContents()
        self.to_de.setEnabled(True)
        self.mw.de_page.counts.setText(r["counts_path"])
        self.mw.de_page.annot.setText(r["annotation_path"])
        self.mw.de_page.load_counts()
        if r.get("excluded"):
            QMessageBox.warning(self, APP_NAME, "Counts completed. BAM files excluded because damaged:\n" +
                                "\n".join(f"• {Path(k).name}: {v}" for k, v in r["excluded"].items()))


# ----------------------------------------------------------------------------- DE
ROLE_LABEL = {"ctrl": "Control", "case": "Case", "other": "In model", "excl": "Excluded"}
ROLE_KEY = {"ctrl": "C", "case": "T", "other": "M", "excl": "E"}
ROLE_FG = {"ctrl": "#1f5a99", "case": "#b02a3f", "other": "#4f5b68", "excl": "#a0a8b3"}
ROLE_BG = {"ctrl": "#e4eef9", "case": "#fbe7ea", "other": None, "excl": "#f3f4f6"}
ROLE_TIP = {"ctrl": "Reference of the comparison (denominator of the fold change)",
            "case": "Group compared with the control (numerator of the fold change)",
            "other": "Not compared, but included in the model: contributes to dispersion estimation",
            "excl": "Excluded from the analysis"}


def task_deseq_batch(task, counts, jobs, cfg):
    results, failed = [], {}
    for i, j in enumerate(jobs):
        task.progress.emit(i, len(jobs))
        task.log.emit(f"━━ Comparison {i + 1}/{len(jobs)}: {j['name']}")
        try:
            r = task_deseq(task, counts, j["groups"], j["batches"], j["ref"], j["case"], cfg)
            r["name"] = j["name"]
            r["job"] = j
            results.append(r)
        except Cancelled:
            raise
        except Exception as e:
            failed[j["name"]] = str(e)
            task.log.emit(f"✖ {j['name']}: {e}")
    task.progress.emit(len(jobs), len(jobs))
    if not results:
        raise RuntimeError("No comparison completed:\n" + "\n".join(f"{k}: {v}" for k, v in failed.items()))
    return {"results": results, "failed": failed}


def make_de_plots(r, top_n):
    res, cfg = r["res"], r["cfg"]
    vst, grp, bat, lfc = r["vst"], r["groups"], r["batches"], r["lfc_col"]
    ref, case, a, t = r["ref"], r["case"], cfg["alpha"], cfg["lfc"]
    return {
        "PCA": lambda f: plot_pca(f, vst, grp, bat, ref, case),
        "Volcano plot": lambda f: plot_volcano(f, res, a, t, lfc, f"{case} vs {ref}"),
        "MA plot": lambda f: plot_ma(f, res, a, t, lfc),
        f"Heatmap top {top_n} genes": lambda f: plot_heatmap(f, vst, res, grp, ref, case, top_n, a),
        "Sample correlation": lambda f: plot_corr(f, vst, grp, ref, case),
    }


def export_plots(plots, outdir):
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    n = 0
    for name, fn in plots.items():
        f = Figure(figsize=(8, 6.5), dpi=100, layout="constrained", facecolor="white")
        try:
            fn(f)
        except Exception:
            continue
        base = outdir / sanitize(name)
        f.savefig(f"{base}.png", dpi=300)
        f.savefig(f"{base}.pdf")
        n += 1
    return n


class DEPage(Page):
    C_NAME, C_ROLE, C_GROUP, C_BATCH, C_READS = range(5)

    def __init__(self, mw):
        super().__init__(mw, "3  Differential expression",
                         "Load a raw count matrix, give each sample a role (control, case, "
                         "in model, excluded), individually or by group, and run one or more DESeq2 comparisons.")
        self.counts_df, self.symbols = None, None
        self.result, self.results, self.queue, self.roles = None, [], [], {}
        self._auto_names = True
        split = QSplitter(Qt.Orientation.Horizontal)
        left = QWidget()
        L = QVBoxLayout(left)
        L.setContentsMargins(0, 0, 0, 0)

        # --- data
        data = card("Data")
        f = QFormLayout(data)
        self.counts = PathEdit(filt="Tables (*.tsv *.txt *.csv *.tsv.gz *.txt.gz *.csv.gz);;All (*)",
                               placeholder="Raw count matrix")
        lb = QPushButton("Load counts")
        lb.clicked.connect(self.load_counts)
        self.annot = PathEdit(filt="GTF or table (*.gtf *.gtf.gz *.tsv *.csv *.txt);;All (*)",
                              placeholder="Optional: GTF or gene_id → symbol table")
        self.meta = PathEdit(filt="Tables (*.tsv *.csv *.txt);;All (*)",
                             placeholder="Optional: sample, group, batch")
        mb = QPushButton("Apply")
        mb.clicked.connect(self.apply_metadata)
        r1, r2 = QHBoxLayout(), QHBoxLayout()
        r1.addWidget(self.counts, 1)
        r1.addWidget(lb)
        r2.addWidget(self.meta, 1)
        r2.addWidget(mb)
        f.addRow("Counts", r1)
        f.addRow("Annotation", self.annot)
        f.addRow("Metadata", r2)
        L.addWidget(data)

        # --- experimental design
        des = card("Samples: controls and cases")
        dl = QVBoxLayout(des)
        top = QHBoxLayout()
        self.search = QLineEdit()
        self.search.setPlaceholderText("Filter samples by name, group or batch (regex allowed)…")
        self.search.textChanged.connect(self.apply_filter)
        selv = QPushButton("Select visible")
        selv.clicked.connect(self.select_visible)
        top.addWidget(self.search, 1)
        top.addWidget(selv)
        dl.addLayout(top)
        dl.addWidget(hint("Groups: click to select their samples (Ctrl+click to add them to the selection)."))
        chipw = QWidget()
        self.chip_box = QHBoxLayout(chipw)
        self.chip_box.setContentsMargins(0, 0, 0, 0)
        self.chip_box.addStretch()
        chips = QScrollArea()
        chips.setWidget(chipw)
        chips.setWidgetResizable(True)
        chips.setFixedHeight(44)
        chips.setVerticalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAlwaysOff)
        chips.setFrameShape(QFrame.Shape.NoFrame)
        dl.addWidget(chips)
        self.table = QTableWidget(0, 5)
        self.table.setHorizontalHeaderLabels(["Sample", "Role", "Group", "Batch", "Reads"])
        hh = self.table.horizontalHeader()
        hh.setSectionResizeMode(self.C_NAME, QHeaderView.ResizeMode.Stretch)
        self.table.verticalHeader().setVisible(False)
        self.table.setSelectionBehavior(QAbstractItemView.SelectionBehavior.SelectRows)
        self.table.setSelectionMode(QAbstractItemView.SelectionMode.ExtendedSelection)
        self.table.setEditTriggers(QAbstractItemView.EditTrigger.DoubleClicked |
                                   QAbstractItemView.EditTrigger.EditKeyPressed)
        self.table.setMinimumHeight(240)
        self.table.setSortingEnabled(False)
        self.table.setContextMenuPolicy(Qt.ContextMenuPolicy.CustomContextMenu)
        self.table.customContextMenuRequested.connect(self.context_menu)
        self.table.itemChanged.connect(self._item_changed)
        self.table.cellDoubleClicked.connect(self._cycle_role)
        dl.addWidget(self.table, 1)
        rb = QHBoxLayout()
        for role in ("ctrl", "case", "other", "excl"):
            b = QPushButton(f"{ROLE_LABEL[role]}  ({ROLE_KEY[role]})")
            b.setObjectName("role_" + role)
            b.setToolTip(ROLE_TIP[role] + f". Shortcut: key {ROLE_KEY[role]} on the table.")
            b.clicked.connect(lambda _=False, r=role: self.set_role(r))
            rb.addWidget(b)
            sc = QShortcut(QKeySequence(ROLE_KEY[role]), self.table)
            sc.setContext(Qt.ShortcutContext.WidgetShortcut)
            sc.activated.connect(lambda r=role: self.set_role(r))
        dl.addLayout(rb)
        tb = QHBoxLayout()
        for txt, fn in (("Assign group…", self.assign_group), ("Set batch…", self.assign_batch),
                        ("Groups from names", self.guess_groups), ("Reset roles", self.reset_roles)):
            b = QPushButton(txt)
            b.clicked.connect(fn)
            tb.addWidget(b)
        dl.addLayout(tb)
        self.design_summary = QLabel("Load a count matrix to start.")
        self.design_summary.setObjectName("designSummary")
        self.design_summary.setWordWrap(True)
        dl.addWidget(self.design_summary)
        L.addWidget(des, 1)

        # --- comparison and queue
        con = card("Comparison")
        cf = QFormLayout(con)
        self.ref_name, self.case_name = QLineEdit(), QLineEdit()
        self.ref_name.setPlaceholderText("e.g. WT")
        self.case_name.setPlaceholderText("e.g. D801N")
        for le in (self.ref_name, self.case_name):
            le.textEdited.connect(lambda _: setattr(self, "_auto_names", False))
            le.textChanged.connect(lambda _: self.update_summary())
        self.only_contrast = QCheckBox("Estimate dispersion from controls and cases only (ignore 'In model' samples)")
        self.use_batch = QCheckBox("Correct for batch (~ batch + condition)")
        cf.addRow("Control name", self.ref_name)
        cf.addRow("Case name", self.case_name)
        cf.addRow(self.only_contrast)
        cf.addRow(self.use_batch)
        qb = QHBoxLayout()
        addq = QPushButton("Add to queue")
        addq.clicked.connect(self.add_to_queue)
        rmq = QPushButton("Remove from queue")
        rmq.clicked.connect(self.remove_from_queue)
        qb.addWidget(addq)
        qb.addWidget(rmq)
        qb.addStretch()
        cf.addRow(qb)
        self.queue_list = QListWidget()
        self.queue_list.setMaximumHeight(96)
        self.queue_list.itemDoubleClicked.connect(self.load_job)
        cf.addRow(self.queue_list)
        cf.addRow(hint("Each queued comparison keeps its own sample selection; double-click to "
                       "reload and edit it."))
        L.addWidget(con)

        # --- parameters
        par = card("Parameters")
        pf = QFormLayout(par)
        self.alpha = QDoubleSpinBox()
        self.alpha.setRange(0.0001, 0.5)
        self.alpha.setDecimals(4)
        self.alpha.setSingleStep(0.01)
        self.alpha.setValue(0.05)
        self.lfc = QDoubleSpinBox()
        self.lfc.setRange(0, 10)
        self.lfc.setSingleStep(0.25)
        self.lfc.setValue(1.0)
        self.min_count = QSpinBox()
        self.min_count.setRange(0, 10000)
        self.min_count.setValue(10)
        self.min_samples = QSpinBox()
        self.min_samples.setRange(0, 1000)
        self.min_samples.setSpecialValueText("auto (smallest group)")
        self.shrink = QCheckBox("LFC shrinkage (apeGLM) for volcano plot and ranking")
        self.shrink.setChecked(True)
        self.top_n = QSpinBox()
        self.top_n.setRange(10, 200)
        self.top_n.setValue(50)
        pf.addRow("padj (α)", self.alpha)
        pf.addRow("Minimum |log2FC|", self.lfc)
        pf.addRow("Minimum reads per gene", self.min_count)
        pf.addRow("…in at least N samples", self.min_samples)
        pf.addRow(self.shrink)
        pf.addRow("Genes in heatmap", self.top_n)
        L.addWidget(par)
        runrow = QHBoxLayout()
        self.run_btn = QPushButton("Run this comparison")
        self.run_btn.setObjectName("primary")
        self.run_btn.clicked.connect(self.run_current)
        self.run_q_btn = QPushButton("Run queue (0)")
        self.run_q_btn.clicked.connect(self.run_queue)
        runrow.addWidget(self.run_btn, 1)
        runrow.addWidget(self.run_q_btn)
        L.addLayout(runrow)
        split.addWidget(scroll(left))

        # --- results
        rightw = QWidget()
        rl = QVBoxLayout(rightw)
        rl.setContentsMargins(0, 0, 0, 0)
        sel = QHBoxLayout()
        sel.addWidget(QLabel("Comparison"))
        self.res_combo = QComboBox()
        self.res_combo.currentIndexChanged.connect(self._select_result)
        sel.addWidget(self.res_combo, 1)
        rl.addLayout(sel)
        right = QTabWidget()
        rt = QWidget()
        rtl = QVBoxLayout(rt)
        self.summary = QLabel("No results.")
        self.summary.setObjectName("summary")
        self.summary.setWordWrap(True)
        rtl.addWidget(self.summary)
        fr = QHBoxLayout()
        self.filter = QLineEdit()
        self.filter.setPlaceholderText("Search gene…")
        self.only_sig = QCheckBox("Significant only")
        self.filter.textChanged.connect(lambda t: self.res_proxy.setFilterFixedString(t))
        self.only_sig.toggled.connect(self.update_table)
        fr.addWidget(self.filter, 1)
        fr.addWidget(self.only_sig)
        rtl.addLayout(fr)
        self.res_view, self.res_model, self.res_proxy = make_table()
        rtl.addWidget(self.res_view, 1)
        br = QHBoxLayout()
        self.open_btn = QPushButton("Open results folder")
        self.open_btn.clicked.connect(lambda: self.result and self.mw.open_path(self.result["outdir"]))
        self.to_pw = QPushButton("Analyse pathways for this comparison")
        self.to_pw.setObjectName("primary")
        self.to_pw.clicked.connect(self.send_to_pathway)
        for b in (self.open_btn, self.to_pw):
            b.setEnabled(False)
            br.addWidget(b)
        br.addStretch()
        rtl.addLayout(br)
        right.addTab(rt, "Results")
        self.plots = PlotPanel()
        right.addTab(self.plots, "Plots")
        self.right_tabs = right
        rl.addWidget(right, 1)
        split.addWidget(rightw)
        split.setSizes([780, 560])
        self.body.addWidget(split, 1)

    # ------------------------------------------------------------ sample table
    def _rows(self):
        return range(self.table.rowCount())

    def _name(self, r):
        return self.table.item(r, self.C_NAME).text()

    def _group(self, r):
        return self.table.item(r, self.C_GROUP).text().strip()

    def selected_rows(self):
        return sorted({i.row() for i in self.table.selectedIndexes()})

    def load_counts(self):
        p = self.counts.text()
        if not p or not Path(p).is_file():
            QMessageBox.warning(self, APP_NAME, "Select a valid count file.")
            return
        try:
            counts, _, symbols, nonint = read_counts(p)
        except Exception as e:
            QMessageBox.critical(self, APP_NAME, f"Cannot read the counts:\n{e}")
            return
        self.counts_df, self.symbols = counts, symbols
        if nonint:
            QMessageBox.warning(self, APP_NAME, "The file contains non-integer values (e.g. Salmon/RSEM estimates): "
                                                "they have been rounded. Do not use TPM/FPKM with DESeq2.")
        self.table.blockSignals(True)
        self.table.setRowCount(0)
        for s in counts.columns:
            r = self.table.rowCount()
            self.table.insertRow(r)
            name = QTableWidgetItem(s)
            name.setFlags(Qt.ItemFlag.ItemIsEnabled | Qt.ItemFlag.ItemIsSelectable)
            role = QTableWidgetItem()
            role.setFlags(Qt.ItemFlag.ItemIsEnabled | Qt.ItemFlag.ItemIsSelectable)
            rd = QTableWidgetItem(f"{int(counts[s].sum()):,}".replace(",", "."))
            rd.setFlags(Qt.ItemFlag.ItemIsEnabled | Qt.ItemFlag.ItemIsSelectable)
            rd.setTextAlignment(Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
            for c, it in ((self.C_NAME, name), (self.C_ROLE, role), (self.C_GROUP, QTableWidgetItem(guess_group(s))),
                          (self.C_BATCH, QTableWidgetItem("")), (self.C_READS, rd)):
                self.table.setItem(r, c, it)
        self.table.blockSignals(False)
        self.roles = {s: "other" for s in counts.columns}
        self._auto_names = True
        self.guess_roles()
        self.table.resizeColumnsToContents()
        self.table.horizontalHeader().setSectionResizeMode(self.C_NAME, QHeaderView.ResizeMode.Stretch)
        self.mw.log(f"Counts loaded: {counts.shape[0]} genes × {counts.shape[1]} samples from {p}")

    def guess_roles(self):
        """Initial proposal: group whose name looks like a control → Control; if there are two groups,
        the other becomes Case. Everything else stays 'In model'."""
        groups = sorted({self._group(r) for r in self._rows()})
        ctrl = next((g for g in groups if CONTROL_RE.search(g)), None)
        others = [g for g in groups if g != ctrl]
        for r in self._rows():
            g, s = self._group(r), self._name(r)
            self.roles[s] = "ctrl" if g == ctrl else ("case" if ctrl and len(others) == 1 else "other")
        self.refresh()

    def refresh(self):
        self.table.blockSignals(True)
        for r in self._rows():
            role = self.roles.get(self._name(r), "other")
            it = self.table.item(r, self.C_ROLE)
            it.setText(ROLE_LABEL[role])
            it.setToolTip(ROLE_TIP[role])
            fg = QBrush(QColor(ROLE_FG[role]))
            bg = QBrush(QColor(ROLE_BG[role])) if ROLE_BG[role] else QBrush()
            for c in range(self.table.columnCount()):
                cell = self.table.item(r, c)
                cell.setBackground(bg)
                cell.setForeground(fg if c in (self.C_ROLE,) or role == "excl" else QBrush(QColor("#1e2a36")))
            f = it.font()
            f.setBold(role in ("ctrl", "case"))
            it.setFont(f)
        self.table.blockSignals(False)
        self.rebuild_chips()
        self._autoname()
        self.update_summary()

    def rebuild_chips(self):
        while self.chip_box.count() > 1:
            w = self.chip_box.takeAt(0).widget()
            if w:
                w.deleteLater()
        counts = {}
        for r in self._rows():
            counts.setdefault(self._group(r) or "(no group)", []).append(r)
        for g in sorted(counts):
            rows = counts[g]
            roles = {self.roles.get(self._name(r)) for r in rows}
            b = QToolButton()
            b.setText(f"{g}  ·  {len(rows)}")
            b.setObjectName("chip_" + (roles.pop() if len(roles) == 1 else "mixed"))
            b.setToolTip("Click: select the samples of this group")
            b.clicked.connect(lambda _=False, rr=rows: self.select_rows(rr))
            self.chip_box.insertWidget(self.chip_box.count() - 1, b)

    def select_rows(self, rows):
        add = bool(QApplication.keyboardModifiers() & Qt.KeyboardModifier.ControlModifier)
        sm = self.table.selectionModel()
        if not add:
            self.table.clearSelection()
        from PySide6.QtCore import QItemSelectionModel
        for r in rows:
            sm.select(self.table.model().index(r, 0),
                      QItemSelectionModel.SelectionFlag.Select | QItemSelectionModel.SelectionFlag.Rows)
        self.table.setFocus()

    def apply_filter(self, text):
        try:
            rx = re.compile(text, re.I) if text else None
        except re.error:
            rx = re.compile(re.escape(text), re.I)
        for r in self._rows():
            hay = " ".join(self.table.item(r, c).text() for c in (self.C_NAME, self.C_GROUP, self.C_BATCH))
            self.table.setRowHidden(r, bool(rx) and not rx.search(hay))

    def select_visible(self):
        self.select_rows([r for r in self._rows() if not self.table.isRowHidden(r)])

    def set_role(self, role, rows=None):
        rows = self.selected_rows() if rows is None else rows
        if not rows:
            self.mw.status.setText("Select one or more samples (also via groups or the filter).")
            return
        for r in rows:
            self.roles[self._name(r)] = role
        self.refresh()

    def _cycle_role(self, row, col):
        if col == self.C_ROLE:
            order = ["ctrl", "case", "other", "excl"]
            cur = self.roles.get(self._name(row), "other")
            self.set_role(order[(order.index(cur) + 1) % 4], [row])

    def _item_changed(self, item):
        if item.column() in (self.C_GROUP, self.C_BATCH):
            self.rebuild_chips()
            self._autoname()
            self.update_summary()

    def context_menu(self, pos):
        rows = self.selected_rows()
        if not rows:
            return
        m = QMenu(self)
        for role in ("ctrl", "case", "other", "excl"):
            m.addAction(f"Set as {ROLE_LABEL[role]}  ({ROLE_KEY[role]})", lambda r=role: self.set_role(r))
        m.addSeparator()
        m.addAction("Select whole group", lambda: self.select_rows(
            [r for r in self._rows() if self._group(r) in {self._group(x) for x in rows}]))
        m.addAction("Assign group…", self.assign_group)
        m.addAction("Set batch…", self.assign_batch)
        m.exec(self.table.viewport().mapToGlobal(pos))

    def _groups_existing(self):
        return sorted({self._group(r) for r in self._rows() if self._group(r)})

    def assign_group(self):
        rows = self.selected_rows()
        if not rows:
            return
        g, ok = QInputDialog.getItem(self, "Group", "Group name:", self._groups_existing() or ["Group1"], 0, True)
        if ok and g.strip():
            for r in rows:
                self.table.item(r, self.C_GROUP).setText(g.strip())

    def assign_batch(self):
        rows = self.selected_rows()
        if not rows:
            return
        b, ok = QInputDialog.getText(self, "Batch", "Batch (e.g. lot, run, extraction date):")
        if ok:
            for r in rows:
                self.table.item(r, self.C_BATCH).setText(b.strip())
            if b.strip():
                self.use_batch.setChecked(True)

    def guess_groups(self):
        self.table.blockSignals(True)
        for r in self._rows():
            self.table.item(r, self.C_GROUP).setText(guess_group(self._name(r)))
        self.table.blockSignals(False)
        self.refresh()

    def reset_roles(self):
        self.roles = {self._name(r): "other" for r in self._rows()}
        self._auto_names = True
        self.ref_name.clear()
        self.case_name.clear()
        self.refresh()

    def _autoname(self):
        if not self._auto_names:
            return
        for role, le, default in (("ctrl", self.ref_name, "Control"), ("case", self.case_name, "Case")):
            gs = {self._group(r) for r in self._rows() if self.roles.get(self._name(r)) == role}
            gs.discard("")
            le.blockSignals(True)
            le.setText(gs.pop() if len(gs) == 1 else (default if gs else ""))
            le.blockSignals(False)

    def update_summary(self):
        if not self.table.rowCount():
            return
        n = {k: sum(1 for v in self.roles.values() if v == k) for k in ROLE_LABEL}
        ref, case = self.ref_name.text().strip() or "control", self.case_name.text().strip() or "case"
        parts = [f"<b style='color:{ROLE_FG['case']}'>{case}</b> ({n['case']}) vs "
                 f"<b style='color:{ROLE_FG['ctrl']}'>{ref}</b> ({n['ctrl']})",
                 f"{n['other']} in model" + (" (ignored)" if self.only_contrast.isChecked() else ""),
                 f"{n['excl']} excluded"]
        notes = []
        if not n["ctrl"] or not n["case"]:
            notes.append("Select at least one control and one case.")
        elif min(n["ctrl"], n["case"]) < 2:
            notes.append("Fewer than 2 replicates per group: the analysis will run, with exploratory results.")
        if ref == case:
            notes.append("Control and case must have different names.")
        html = " · ".join(parts)
        if notes:
            html += "<br><span style='color:#b4690e'>" + " ".join(notes) + "</span>"
        self.design_summary.setText(html)

    def apply_metadata(self):
        p = self.meta.text()
        if not p or self.table.rowCount() == 0:
            QMessageBox.information(self, APP_NAME, "Load the counts first and select the metadata file.")
            return
        try:
            md = pd.read_csv(p, sep=_sniff_sep(p), dtype=str)
        except Exception as e:
            QMessageBox.critical(self, APP_NAME, f"Unreadable metadata:\n{e}")
            return
        cols = {c.lower(): c for c in md.columns}
        sc = next((cols[c] for c in ("sample", "sample", "sample_id", "id", "name") if c in cols), md.columns[0])
        gc = next((cols[c] for c in ("condition", "group", "group", "condizione", "treatment", "trattamento")
                   if c in cols), md.columns[1] if len(md.columns) > 1 else None)
        bc = next((cols[c] for c in ("batch", "lotto", "run") if c in cols), None)
        rc = next((cols[c] for c in ("role", "ruolo") if c in cols), None)
        md[sc] = md[sc].map(clean_sample_name)
        md = md.set_index(sc)
        role_map = {"controllo": "ctrl", "control": "ctrl", "ctrl": "ctrl", "caso": "case", "case": "case",
                    "escluso": "excl", "exclude": "excl", "excluded": "excl", "excl": "excl"}
        hit = 0
        self.table.blockSignals(True)
        for r in self._rows():
            s = self._name(r)
            if s in md.index:
                hit += 1
                if gc:
                    self.table.item(r, self.C_GROUP).setText(str(md.at[s, gc]))
                if bc:
                    self.table.item(r, self.C_BATCH).setText(str(md.at[s, bc]))
                if rc:
                    self.roles[s] = role_map.get(str(md.at[s, rc]).strip().lower(), "other")
            else:
                self.roles[s] = "excl"
        self.table.blockSignals(False)
        if not rc:
            self.guess_roles()
            for r in self._rows():
                if self._name(r) not in md.index:
                    self.roles[self._name(r)] = "excl"
        if bc:
            self.use_batch.setChecked(True)
        self.refresh()
        QMessageBox.information(self, APP_NAME, f"Metadata applied to {hit} samples out of {self.table.rowCount()}; "
                                                "samples missing from the metadata have been excluded.")

    # ------------------------------------------------------------ comparisons
    def collect_job(self, quiet=False):
        ref, case = self.ref_name.text().strip(), self.case_name.text().strip()
        ctrl = [s for s, r in self.roles.items() if r == "ctrl"]
        cas = [s for s, r in self.roles.items() if r == "case"]
        err = None
        if not ctrl or not cas:
            err = "Assign at least one sample as Control and one as Case."
        elif not ref or not case or ref == case:
            err = "Provide two different names for control and case."
        if err:
            if not quiet:
                QMessageBox.warning(self, APP_NAME, err)
            return None
        groups, batches = {}, {}
        for r in self._rows():
            s, role = self._name(r), self.roles.get(self._name(r))
            if role == "ctrl":
                groups[s] = ref
            elif role == "case":
                groups[s] = case
            elif role == "other" and not self.only_contrast.isChecked():
                g = self._group(r) or "others"
                groups[s] = g if g not in (ref, case) else f"{g}_others"
            else:
                continue
            batches[s] = self.table.item(r, self.C_BATCH).text().strip()
        return {"name": f"{case} vs {ref}", "ref": ref, "case": case, "groups": groups, "batches": batches,
                "roles": dict(self.roles), "only_contrast": self.only_contrast.isChecked(),
                "use_batch": self.use_batch.isChecked(),
                "group_col": {self._name(r): self._group(r) for r in self._rows()}}

    def add_to_queue(self):
        job = self.collect_job()
        if not job:
            return
        self.queue = [j for j in self.queue if j["name"] != job["name"]] + [job]
        self._refresh_queue()

    def remove_from_queue(self):
        names = {i.data(Qt.ItemDataRole.UserRole) for i in self.queue_list.selectedItems()}
        self.queue = [j for j in self.queue if j["name"] not in names]
        self._refresh_queue()

    def _refresh_queue(self):
        self.queue_list.clear()
        for j in self.queue:
            nc = sum(1 for v in j["roles"].values() if v == "ctrl")
            nk = sum(1 for v in j["roles"].values() if v == "case")
            it = QListWidgetItem(f"{j['name']}   ({nk} vs {nc}; {len(j['groups'])} samples in model)")
            it.setData(Qt.ItemDataRole.UserRole, j["name"])
            self.queue_list.addItem(it)
        self.run_q_btn.setText(f"Run queue ({len(self.queue)})")

    def load_job(self, item):
        j = next((x for x in self.queue if x["name"] == item.data(Qt.ItemDataRole.UserRole)), None)
        if not j:
            return
        self.table.blockSignals(True)
        for r in self._rows():
            s = self._name(r)
            if s in j["group_col"]:
                self.table.item(r, self.C_GROUP).setText(j["group_col"][s])
            if s in j["batches"]:
                self.table.item(r, self.C_BATCH).setText(j["batches"][s])
        self.table.blockSignals(False)
        self.roles = {s: j["roles"].get(s, "other") for s in self.roles}
        self._auto_names = False
        self.ref_name.setText(j["ref"])
        self.case_name.setText(j["case"])
        self.only_contrast.setChecked(j["only_contrast"])
        self.use_batch.setChecked(j["use_batch"])
        self.refresh()

    def _cfg(self, use_batch):
        return dict(outdir=self.mw.outdir(), threads=self.mw.threads(), alpha=self.alpha.value(), lfc=self.lfc.value(),
                    min_count=self.min_count.value(), min_samples=self.min_samples.value(),
                    shrink=self.shrink.isChecked(), only_contrast=False, use_batch=use_batch,
                    annotation=self.annot.text(), symbols=self.symbols, counts_file=self.counts.text())

    def _start(self, jobs):
        if self.counts_df is None:
            QMessageBox.warning(self, APP_NAME, "Load the count matrix first.")
            return
        use_batch = any(j["use_batch"] for j in jobs)
        title = jobs[0]["name"] if len(jobs) == 1 else f"{len(jobs)} DESeq2 comparisons"
        self.mw.start_task(task_deseq_batch, self.counts_df, jobs, self._cfg(use_batch), title=title,
                           on_done=self.done)

    def run_current(self):
        job = self.collect_job()
        if job:
            self._start([job])

    def run_queue(self):
        if not self.queue:
            job = self.collect_job()
            if job:
                self._start([job])
            return
        self._start(list(self.queue))

    def done(self, out):
        warns, top_n = [], self.top_n.value()
        for r in out["results"]:
            self.results = [x for x in self.results if x["name"] != r["name"]] + [r]
            try:
                export_plots(make_de_plots(r, top_n), Path(r["outdir"]) / "plots")
            except Exception as e:
                self.mw.log(f"Plot export failed for {r['name']}: {e}")
            warns += [f"{r['name']}: {w}" for w in r.get("warnings", [])]
        self.res_combo.blockSignals(True)
        self.res_combo.clear()
        self.res_combo.addItems([x["name"] for x in self.results])
        self.res_combo.blockSignals(False)
        self.res_combo.setCurrentIndex(len(self.results) - 1)
        self._select_result(len(self.results) - 1)
        msgs = []
        if out["failed"]:
            msgs.append("Comparisons not completed:\n" + "\n".join(f"• {k}: {v}" for k, v in out["failed"].items()))
        if warns:
            msgs.append("Warnings:\n" + "\n".join(f"• {w}" for w in warns))
        if msgs:
            QMessageBox.warning(self, APP_NAME, "\n\n".join(msgs))

    def _select_result(self, idx):
        if idx < 0 or idx >= len(self.results):
            return
        r = self.result = self.results[idx]
        res, cfg = r["res"], r["cfg"]
        ngenes = f"{len(res):,}".replace(",", ".")
        txt = (f"{r['case']} vs {r['ref']}: {ngenes} genes tested — {r['n_up']} UP, {r['n_down']} DOWN "
               f"(padj < {cfg['alpha']}, |log2FC| ≥ {cfg['lfc']}).")
        if r.get("mode") == "exploratory":
            txt += " Exploratory analysis: threshold on log2FC only."
        if r.get("warnings"):
            txt = "⚠ " + txt
        self.summary.setText(txt)
        self.update_table()
        self.plots.set_plots(make_de_plots(r, self.top_n.value()), str(Path(r["outdir"]) / "plots"))
        for b in (self.open_btn, self.to_pw):
            b.setEnabled(True)

    def update_table(self):
        if not self.result:
            return
        res = self.result["res"]
        if self.only_sig.isChecked():
            res = res[res["direction"] != ""]
        self.res_model.set_df(res.reset_index())
        self.res_view.resizeColumnsToContents()

    def send_to_pathway(self):
        if self.result:
            self.mw.pw_page.set_de(self.result)
            self.mw.goto(self.mw.pw_page)


# ----------------------------------------------------------------------------- Pathway
class EnrichmentDialog(QDialog):
    def __init__(self, parent, ranking, info, term, lib):
        super().__init__(parent)
        self.setWindowTitle(f"GSEA — {term}")
        self.resize(820, 620)
        lay = QVBoxLayout(self)
        fig = Figure(figsize=(8, 6), dpi=100, facecolor="white")
        canvas = FigureCanvas(fig)
        lay.addWidget(NavToolbar(canvas, self))
        lay.addWidget(canvas, 1)
        plot_enrichment(fig, ranking, info, term, lib)
        canvas.draw()


class PathwayPage(Page):
    def __init__(self, mw):
        super().__init__(mw, "4  Pathway", "Over-representation analysis (hypergeometric test, universe = tested genes) "
                                          "on UP and DOWN genes, and preranked GSEA on the full ranked list. "
                                          "Enrichr libraries (KEGG, Reactome, GO, Hallmark, WikiPathways) or custom "
                                          "GMT files, e.g. from MSigDB.")
        self.de, self.result = None, None
        split = QSplitter(Qt.Orientation.Horizontal)
        left = QWidget()
        L = QVBoxLayout(left)
        L.setContentsMargins(0, 0, 0, 0)

        src = card("Input")
        sf = QVBoxLayout(src)
        self.src_session = QRadioButton("DE results from this session")
        self.src_file = QRadioButton("DE results file (TSV/CSV with log2FoldChange and padj)")
        self.src_file.setChecked(True)
        self.src_label = hint("No DE results in memory.")
        self.de_file = PathEdit(filt="Tables (*.tsv *.csv *.txt);;All (*)", placeholder="DE_results.tsv")
        sf.addWidget(self.src_session)
        sf.addWidget(self.src_label)
        sf.addWidget(self.src_file)
        sf.addWidget(self.de_file)
        L.addWidget(src)

        gsb = card("Gene set")
        gl = QVBoxLayout(gsb)
        orow = QHBoxLayout()
        self.organism = QComboBox()
        self.organism.addItems(["Human", "Mouse"])
        self.organism.currentTextChanged.connect(self.reset_libs)
        upd = QPushButton("Update list from Enrichr")
        upd.clicked.connect(self.fetch_libs)
        orow.addWidget(QLabel("Organism"))
        orow.addWidget(self.organism)
        orow.addStretch()
        orow.addWidget(upd)
        gl.addLayout(orow)
        self.lib_filter = QLineEdit()
        self.lib_filter.setPlaceholderText("Filter libraries…")
        self.lib_filter.textChanged.connect(self.filter_libs)
        gl.addWidget(self.lib_filter)
        self.libs = QListWidget()
        self.libs.setMinimumHeight(170)
        gl.addWidget(self.libs, 1)
        grow = QHBoxLayout()
        addg = QPushButton("Add GMT file…")
        addg.clicked.connect(self.add_gmt)
        clg = QPushButton("Remove GMT")
        clg.clicked.connect(lambda: [self.gmts.takeItem(self.gmts.row(i)) for i in self.gmts.selectedItems()])
        grow.addWidget(addg)
        grow.addWidget(clg)
        grow.addStretch()
        gl.addLayout(grow)
        self.gmts = QListWidget()
        self.gmts.setMaximumHeight(70)
        gl.addWidget(self.gmts)
        gl.addWidget(hint(f"Downloaded libraries are cached in {GENESET_CACHE} and can then be used offline."))
        L.addWidget(gsb, 1)

        mt = card("Methods")
        mf = QFormLayout(mt)
        self.do_ora = QCheckBox("ORA on UP and DOWN genes")
        self.do_ora.setChecked(True)
        self.do_gsea = QCheckBox("GSEA preranked")
        self.do_gsea.setChecked(True)
        self.alpha = QDoubleSpinBox()
        self.alpha.setRange(0.0001, 0.5)
        self.alpha.setDecimals(4)
        self.alpha.setValue(0.05)
        self.lfc = QDoubleSpinBox()
        self.lfc.setRange(0, 10)
        self.lfc.setSingleStep(0.25)
        self.lfc.setValue(1.0)
        self.metric = QComboBox()
        self.metric.addItems(["Wald statistic (stat)", "sign(log2FC) × −log10(p)", "log2FC (shrunk if available)"])
        self.perm = QSpinBox()
        self.perm.setRange(100, 10000)
        self.perm.setSingleStep(100)
        self.perm.setValue(1000)
        self.min_size = QSpinBox()
        self.min_size.setRange(3, 500)
        self.min_size.setValue(15)
        self.max_size = QSpinBox()
        self.max_size.setRange(20, 5000)
        self.max_size.setValue(500)
        mf.addRow(self.do_ora)
        mf.addRow("ORA: padj <", self.alpha)
        mf.addRow("ORA: |log2FC| ≥", self.lfc)
        mf.addRow(self.do_gsea)
        mf.addRow("GSEA: metric", self.metric)
        mf.addRow("GSEA: permutations", self.perm)
        mf.addRow("Minimum gene-set size", self.min_size)
        mf.addRow("Maximum gene-set size", self.max_size)
        L.addWidget(mt)
        run = QPushButton("Run pathway analysis")
        run.setObjectName("primary")
        run.clicked.connect(self.run)
        L.addWidget(run)
        split.addWidget(scroll(left))

        right = QTabWidget()
        ow = QWidget()
        ol = QVBoxLayout(ow)
        self.ora_filter = QLineEdit()
        self.ora_filter.setPlaceholderText("Search term or gene…")
        ol.addWidget(self.ora_filter)
        self.ora_view, self.ora_model, self.ora_proxy = make_table()
        self.ora_filter.textChanged.connect(self.ora_proxy.setFilterFixedString)
        ol.addWidget(self.ora_view, 1)
        right.addTab(ow, "ORA")
        gw = QWidget()
        glw = QVBoxLayout(gw)
        self.gsea_filter = QLineEdit()
        self.gsea_filter.setPlaceholderText("Search term…")
        glw.addWidget(self.gsea_filter)
        self.gsea_view, self.gsea_model, self.gsea_proxy = make_table()
        self.gsea_filter.textChanged.connect(self.gsea_proxy.setFilterFixedString)
        self.gsea_view.doubleClicked.connect(self.show_enrichment)
        glw.addWidget(self.gsea_view, 1)
        glw.addWidget(hint("Double-click a term for its enrichment plot (running ES)."))
        right.addTab(gw, "GSEA")
        self.plots = PlotPanel()
        right.addTab(self.plots, "Plots")
        split.addWidget(right)
        split.setSizes([560, 760])
        self.body.addWidget(split, 1)
        self.reset_libs("Human")

    def reset_libs(self, organism, names=None):
        checked = {self.libs.item(i).text() for i in range(self.libs.count())
                   if self.libs.item(i).checkState() == Qt.CheckState.Checked} if names else DEFAULT_CHECKED
        self.libs.clear()
        for n in (names or DEFAULT_LIBS.get(organism, [])):
            it = QListWidgetItem(n)
            it.setFlags(it.flags() | Qt.ItemFlag.ItemIsUserCheckable)
            it.setCheckState(Qt.CheckState.Checked if n in checked else Qt.CheckState.Unchecked)
            self.libs.addItem(it)
        self.filter_libs(self.lib_filter.text())

    def filter_libs(self, text):
        t = text.lower()
        for i in range(self.libs.count()):
            it = self.libs.item(i)
            it.setHidden(bool(t) and t not in it.text().lower() and it.checkState() != Qt.CheckState.Checked)

    def fetch_libs(self):
        org = self.organism.currentText()

        def fn(task):
            import gseapy as gp
            task.log.emit(f"Retrieving the list of Enrichr libraries ({org})…")
            return gp.get_library_name(organism=org)
        self.mw.start_task(fn, title="Library list", on_done=lambda names: self.reset_libs(org, sorted(names)))

    def add_gmt(self):
        fs, _ = QFileDialog.getOpenFileNames(self, "GMT file", str(Path.home()), "GMT (*.gmt *.gmt.gz);;All (*)")
        for f in fs:
            if not self.gmts.findItems(f, Qt.MatchFlag.MatchExactly):
                self.gmts.addItem(f)

    def set_de(self, de):
        self.de = de
        self.src_session.setChecked(True)
        self.src_label.setText(f"{de['case']} vs {de['ref']} — {de['n_up']} UP, {de['n_down']} DOWN")
        self.alpha.setValue(de["cfg"]["alpha"])
        self.lfc.setValue(de["cfg"]["lfc"])

    def run(self):
        if self.src_session.isChecked():
            if not self.de:
                QMessageBox.warning(self, APP_NAME, "No DE results in this session: load one from file.")
                return
            res, lfc_col = self.de["res"], self.de["lfc_col"]
            outdir = Path(self.de["outdir"]) / "pathway"
        else:
            p = self.de_file.text()
            if not p or not Path(p).is_file():
                QMessageBox.warning(self, APP_NAME, "Select a DE results file.")
                return
            try:
                res = pd.read_csv(p, sep=_sniff_sep(p), index_col=0)
            except Exception as e:
                QMessageBox.critical(self, APP_NAME, f"Unreadable file:\n{e}")
                return
            needed = {"log2FoldChange", "padj"}
            if not needed.issubset(res.columns):
                QMessageBox.critical(self, APP_NAME, "The file must contain the log2FoldChange and padj columns "
                                                     "(output DESeq2/PyDESeq2).")
                return
            if "symbol" not in res.columns:
                alt = next((c for c in res.columns if str(c).lower() in SYMBOL_COLS), None)
                if alt:
                    res = res.rename(columns={alt: "symbol"})
            for c in ("pvalue", "stat"):
                if c not in res.columns:
                    res[c] = np.nan
            if res["pvalue"].isna().all():
                res["pvalue"] = res["padj"]
            lfc_col = "log2FoldChange_shrunk" if "log2FoldChange_shrunk" in res.columns else "log2FoldChange"
            outdir = Path(self.mw.outdir()) / "05_pathway" / sanitize(Path(p).stem)
        libs = [self.libs.item(i).text() for i in range(self.libs.count())
                if self.libs.item(i).checkState() == Qt.CheckState.Checked]
        gmts = [self.gmts.item(i).text() for i in range(self.gmts.count())]
        if not (self.do_ora.isChecked() or self.do_gsea.isChecked()):
            QMessageBox.warning(self, APP_NAME, "Select at least one method.")
            return
        if not libs and not gmts:
            QMessageBox.warning(self, APP_NAME, "Select at least one library or GMT file.")
            return
        metric = ["stat", "signed_p", "lfc"][self.metric.currentIndex()]
        if metric == "stat" and res["stat"].isna().all():
            metric = "signed_p"
        if metric == "signed_p" and res["pvalue"].isna().all():
            metric = "lfc"
        if self.do_ora.isChecked() and res["padj"].isna().all():
            QMessageBox.information(self, APP_NAME, "The results have no p-values (exploratory analysis): "
                                                    "ORA will find no significant genes; GSEA uses log2FC.")
        cfg = dict(outdir=str(outdir), threads=self.mw.threads(), organism=self.organism.currentText(),
                   libraries=libs, gmts=gmts, ora=self.do_ora.isChecked(), gsea=self.do_gsea.isChecked(),
                   alpha=self.alpha.value(), lfc=self.lfc.value(), lfc_col=lfc_col, metric=metric,
                   perm=self.perm.value(), min_size=self.min_size.value(), max_size=self.max_size.value())
        self.mw.start_task(task_pathway, res, cfg, title="Pathway analysis", on_done=self.done)

    def done(self, r):
        self.result = r
        self.ora_model.set_df(r["ora"])
        self.gsea_model.set_df(r["gsea"])
        for v in (self.ora_view, self.gsea_view):
            v.resizeColumnsToContents()
            for c in range(v.model().columnCount()):
                if v.columnWidth(c) > 380:
                    v.setColumnWidth(c, 380)
        ora_df, gsea_df = r["ora"], r["gsea"]
        plots = {}
        if len(ora_df):
            plots["ORA: top terms"] = lambda f: plot_ora_bars(f, ora_df)
        if len(gsea_df):
            plots["GSEA: NES by direction"] = lambda f: plot_gsea_bars(f, gsea_df)
        plot_dir = Path(r["outdir"]) / "plots"
        self.plots.set_plots(plots, str(plot_dir))
        if plots:
            self.plots.export_all(plot_dir)
        QMessageBox.information(self, APP_NAME, f"Analysis completed. Results in:\n{r['outdir']}")

    def show_enrichment(self, idx):
        if not self.result or self.result["ranking"] is None:
            return
        row = self.gsea_proxy.mapToSource(idx).row()
        df = self.gsea_model._df
        term, lib = df.iloc[row]["Term"], df.iloc[row]["Library"]
        info = self.result["gsea_obj"].get(lib, {}).get(term)
        if not info:
            QMessageBox.information(self, APP_NAME, "Details not available for this term.")
            return
        try:
            EnrichmentDialog(self, self.result["ranking"], info, term, lib).exec()
        except Exception as e:
            QMessageBox.warning(self, APP_NAME, str(e))


# =============================================================================
# 6. Main window
# =============================================================================
class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle(f"{APP_NAME} {APP_VERSION}")
        self.resize(1440, 920)
        self.settings = QSettings(APP_NAME, APP_NAME)
        self.task = None
        self._on_done = None

        central = QWidget()
        self.setCentralWidget(central)
        root = QVBoxLayout(central)
        root.setContentsMargins(0, 0, 0, 0)
        root.setSpacing(0)

        header = QFrame()
        header.setObjectName("header")
        hl = QHBoxLayout(header)
        hl.setContentsMargins(18, 10, 18, 10)
        title = QLabel(APP_NAME)
        title.setObjectName("appTitle")
        sub = QLabel("From reads to pathways, starting from any step")
        sub.setObjectName("appSub")
        self._outdir = PathEdit("dir", placeholder="Project folder")
        self._outdir.setMinimumWidth(420)
        self._outdir.setText(self.settings.value("outdir", str(Path.home() / "RNAseqStudio_project")))
        self._threads = QSpinBox()
        self._threads.setRange(1, os.cpu_count() or 1)
        self._threads.setValue(int(self.settings.value("threads", max(1, (os.cpu_count() or 2) - 2))))
        hl.addWidget(title)
        hl.addSpacing(12)
        hl.addWidget(sub)
        hl.addStretch()
        hl.addWidget(QLabel("Project"))
        hl.addWidget(self._outdir)
        hl.addSpacing(10)
        hl.addWidget(QLabel("Threads"))
        hl.addWidget(self._threads)
        root.addWidget(header)

        vsplit = QSplitter(Qt.Orientation.Vertical)
        top = QWidget()
        tl = QHBoxLayout(top)
        tl.setContentsMargins(0, 0, 0, 0)
        tl.setSpacing(0)
        self.nav = QListWidget()
        self.nav.setObjectName("nav")
        self.nav.setFixedWidth(260)
        self.nav.setHorizontalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAlwaysOff)
        self.stack = QStackedWidget()
        tl.addWidget(self.nav)
        tl.addWidget(self.stack, 1)
        vsplit.addWidget(top)

        logw = QWidget()
        logw.setObjectName("logPanel")
        ll = QVBoxLayout(logw)
        ll.setContentsMargins(14, 8, 14, 8)
        bar = QHBoxLayout()
        self.status = QLabel("Ready")
        self.progress = QProgressBar()
        self.progress.setMaximumWidth(320)
        self.progress.setVisible(False)
        self.cancel_btn = QPushButton("Cancel")
        self.cancel_btn.setEnabled(False)
        self.cancel_btn.clicked.connect(self.cancel_task)
        clr = QPushButton("Clear log")
        bar.addWidget(self.status, 1)
        bar.addWidget(self.progress)
        bar.addWidget(self.cancel_btn)
        bar.addWidget(clr)
        ll.addLayout(bar)
        self.logview = QPlainTextEdit()
        self.logview.setObjectName("log")
        self.logview.setReadOnly(True)
        self.logview.setMaximumBlockCount(20000)
        clr.clicked.connect(self.logview.clear)
        ll.addWidget(self.logview, 1)
        vsplit.addWidget(logw)
        vsplit.setSizes([700, 190])
        root.addWidget(vsplit, 1)

        self.env_page = EnvPage(self)
        self.align_page = AlignPage(self)
        self.count_page = CountPage(self)
        self.de_page = DEPage(self)
        self.pw_page = PathwayPage(self)
        for label, page in (("Environment", self.env_page), ("1  Alignment", self.align_page),
                            ("2  Counts", self.count_page), ("3  Differential expression", self.de_page),
                            ("4  Pathway", self.pw_page)):
            self.nav.addItem(QListWidgetItem(label))
            self.stack.addWidget(page)
        self.nav.currentRowChanged.connect(self.stack.setCurrentIndex)
        missing_tools = [t for t in EXT_TOOLS if not tool_path(t)]
        self.nav.setCurrentRow(3 if missing_tools else 1)
        self.log(f"{APP_NAME} {APP_VERSION} started. Python {platform.python_version()}, "
                 f"PyDESeq2 {_pkg_version('pydeseq2')}, gseapy {_pkg_version('gseapy')}.")
        if missing_tools:
            self.log("Missing external tools (" + ", ".join(missing_tools) + "): alignment and counting "
                     "require installation from the Environment page. DE and pathways are already usable.")

    # --- accessors
    def outdir(self):
        d = self._outdir.text() or str(Path.home() / "RNAseqStudio_project")
        Path(d).mkdir(parents=True, exist_ok=True)
        return d

    def threads(self):
        return self._threads.value()

    def goto(self, page):
        self.nav.setCurrentRow(self.stack.indexOf(page))

    def open_path(self, p):
        QDesktopServices.openUrl(QUrl.fromLocalFile(str(p)))

    def log(self, msg):
        line = f"[{datetime.now():%H:%M:%S}] {msg}"
        self.logview.appendPlainText(line)
        d = self._outdir.text()
        if d and Path(d).is_dir():
            try:
                with open(Path(d) / "rnaseqstudio.log", "a", encoding="utf-8") as fh:
                    fh.write(line + "\n")
            except OSError:
                pass

    # --- task
    def start_task(self, fn, *args, title="", on_done=None):
        if self.task is not None and self.task.isRunning():
            QMessageBox.warning(self, APP_NAME, "Another operation is running.")
            return
        self.task = Task(fn, *args)
        self._on_done, self._title = on_done, title
        self.task.log.connect(self.log)
        self.task.progress.connect(self._progress)
        self.task.done.connect(self._task_done)
        self.task.failed.connect(self._task_failed)
        self.status.setText(f"Running: {title}")
        self.progress.setVisible(True)
        self.progress.setRange(0, 0)
        self.cancel_btn.setEnabled(True)
        self.log(f"▶ {title}")
        self.task.start()

    def _progress(self, i, n):
        self.progress.setRange(0, n)
        self.progress.setValue(i)
        self.progress.setFormat(f"{i}/{n}")

    def _finish(self, ok):
        self.progress.setVisible(False)
        self.cancel_btn.setEnabled(False)
        self.status.setText(("Completed: " if ok else "Interrupted: ") + self._title)

    def _task_done(self, result):
        self._finish(True)
        self.log(f"✔ {self._title} completed")
        if self._on_done:
            try:
                self._on_done(result)
            except Exception as e:
                self.log(traceback.format_exc())
                QMessageBox.critical(self, APP_NAME, f"Error while displaying the results:\n{e}")

    def _task_failed(self, msg):
        self._finish(False)
        self.log(f"✖ {msg}")
        QMessageBox.critical(self, APP_NAME, msg)

    def cancel_task(self):
        if self.task is not None and self.task.isRunning():
            self.task.cancel()
            self.log("Cancellation requested: the running external process is being terminated. "
                     "Python computations (DESeq2, GSEA) stop at the end of the current step.")

    def closeEvent(self, e):
        self.settings.setValue("outdir", self._outdir.text())
        self.settings.setValue("threads", self._threads.value())
        if self.task is not None and self.task.isRunning():
            if QMessageBox.question(self, APP_NAME, "An operation is running. Quit anyway?") != \
                    QMessageBox.StandardButton.Yes:
                e.ignore()
                return
            self.task.cancel()
            self.task.wait(3000)
        e.accept()


QT_LIBS_ENV = APP_HOME / "qtlibs"
# missing library -> conda-forge package that provides it
QT_LIB_PACKAGES = {
    "libxcb-cursor.so.0": "xcb-util-cursor", "libxcb-icccm.so.4": "xcb-util-wm",
    "libxcb-image.so.0": "xcb-util-image", "libxcb-keysyms.so.1": "xcb-util-keysyms",
    "libxcb-render-util.so.0": "xcb-util-renderutil", "libxcb-util.so.1": "xcb-util",
    "libxkbcommon.so.0": "libxkbcommon", "libxkbcommon-x11.so.0": "libxkbcommon",
    "libX11.so.6": "xorg-libx11", "libX11-xcb.so.1": "xorg-libx11", "libfontconfig.so.1": "fontconfig",
    "libdbus-1.so.3": "dbus", "libEGL.so.1": "libegl", "libGL.so.1": "libgl", "libOpenGL.so.0": "libopengl",
}


def _qt_missing_libs(plugin):
    try:
        out = subprocess.run(["ldd", str(plugin)], capture_output=True, text=True).stdout
    except OSError:
        return []
    return sorted({l.split("=>")[0].strip() for l in out.splitlines() if "not found" in l})


def _preload(libdir, names):
    import ctypes
    for n in names:
        f = Path(libdir) / n
        if f.exists():
            ctypes.CDLL(str(f), mode=ctypes.RTLD_GLOBAL)


def prepare_qt_platform():
    """On Linux, check that Qt's xcb platform plugin can be loaded; if system libraries are missing
    (typically libxcb-cursor0), install them from conda-forge into the user's home and preload them,
    without sudo. Also handles the absence of a display (e.g. ssh without -X)."""
    if not sys.platform.startswith("linux") or os.environ.get("QT_QPA_PLATFORM", "").startswith(("offscreen", "vnc")):
        return
    has_x, has_wl = bool(os.environ.get("DISPLAY")), bool(os.environ.get("WAYLAND_DISPLAY"))
    if not has_x and not has_wl:
        sys.exit(f"[{APP_NAME}] No graphical display available (DISPLAY variable is empty).\n"
                 "  - remotely: connect with 'ssh -X user@server' (or -Y), or use X2Go/VNC/NoMachine;\n"
                 "  - on the server itself: start the program from a desktop session.")
    import PySide6
    qt = Path(PySide6.__file__).parent / "Qt"
    plugin = qt / "plugins" / "platforms" / "libqxcb.so"
    if not has_x or not plugin.exists():
        return
    missing = _qt_missing_libs(plugin)
    if not missing:
        return
    _say("Missing system libraries for Qt graphics: " + ", ".join(missing))
    libdir = QT_LIBS_ENV / "lib"
    still = [m for m in missing if not (libdir / m).exists()]
    if still:
        pkgs = sorted({QT_LIB_PACKAGES.get(m, "libxcb" if m.startswith("libxcb") else None) for m in still} - {None})
        unknown = [m for m in still if m not in QT_LIB_PACKAGES and not m.startswith("libxcb")]
        try:
            if not pkgs:
                raise RuntimeError("no known package for " + ", ".join(unknown))
            mm = get_micromamba(_say)
            _say(f"Installing them without sudo from conda-forge into {QT_LIBS_ENV} ({', '.join(pkgs)})…")
            action = "install" if (QT_LIBS_ENV / "conda-meta").exists() else "create"
            subprocess.check_call([mm, action, "-y", "-p", str(QT_LIBS_ENV), "-c", "conda-forge", *pkgs],
                                  env=dict(os.environ, MAMBA_ROOT_PREFIX=str(APP_HOME / "mamba_root")))
        except Exception as e:
            if has_wl:
                _say(f"xcb repair failed ({e}); trying the Wayland backend.")
                os.environ["QT_QPA_PLATFORM"] = "wayland"
                return
            sys.exit(f"[{APP_NAME}] Unable to set up Qt graphics ({e}).\n"
                     f"Solution with privileges: sudo apt install libxcb-cursor0 libxkbcommon-x11-0 "
                     f"libxcb-icccm4 libxcb-image0 libxcb-keysyms1 libxcb-render-util0")
    _preload(libdir, missing)
    _say("Graphics libraries loaded from the application environment.")


def cli():
    """Entry point of the 'rnaseq-studio' command (pip/pipx installation)."""
    if "--version" in sys.argv:
        print(f"{APP_NAME} {APP_VERSION}")
        return
    if "--install-tools" in sys.argv:
        cli_install()
        return
    main()


def main():
    prepare_qt_platform()
    app = QApplication(sys.argv)
    app.setApplicationName(APP_NAME)
    app.setStyle("Fusion")
    apply_palette(app)
    app.setStyleSheet(QSS)
    w = MainWindow()
    w.show()
    sys.exit(app.exec())


if __name__ == "__main__":
    main()
