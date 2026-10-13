"""Unit tests for RNAseqStudio core functions (no external tools or display required)."""
import importlib.util
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

os.environ.setdefault("RNASEQSTUDIO_NO_BOOTSTRAP", "1")
os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("rnaseq_studio", ROOT / "rnaseq_studio.py")
rs = importlib.util.module_from_spec(spec)
sys.modules["rnaseq_studio"] = rs
spec.loader.exec_module(rs)


def test_read_featurecounts_format(tmp_path):
    df = pd.DataFrame({"Chr": "1", "Start": 1, "End": 100, "Strand": "+", "Length": [1000, 2000],
                       "/data/S1.sorted.bam": [10, 0], "/data/S2Aligned.sortedByCoord.out.bam": [5, 7]},
                      index=pd.Index(["ENSG01.3", "ENSG02.1"], name="Geneid"))
    f = tmp_path / "fc.txt"
    with open(f, "w") as fh:
        fh.write("# Program:featureCounts\n")
        df.to_csv(fh, sep="\t")
    counts, lengths, _, nonint = rs.read_counts(f)
    assert list(counts.columns) == ["S1", "S2"]
    assert lengths.tolist() == [1000, 2000]
    assert not nonint and counts.values.sum() == 22


def test_read_counts_rejects_negative(tmp_path):
    f = tmp_path / "c.tsv"
    pd.DataFrame({"A": [1, -2], "B": [3, 4]}, index=["g1", "g2"]).to_csv(f, sep="\t")
    with pytest.raises(ValueError):
        rs.read_counts(f)


def test_discover_fastq(tmp_path):
    for n in ("A_S1_L001_R1_001.fastq.gz", "A_S1_L001_R2_001.fastq.gz", "A_S1_L002_R1_001.fastq.gz",
              "A_S1_L002_R2_001.fastq.gz", "B.fastq"):
        (tmp_path / n).touch()
    found = rs.discover_fastq(tmp_path)
    assert len(found["A"]["R1"]) == 2 and len(found["A"]["R2"]) == 2
    assert found["B"]["R2"] == []


def test_guess_group():
    assert rs.guess_group("WT_1") == "WT"
    assert rs.guess_group("D801N_rep2") == "D801N"


def test_bh_monotone():
    p = np.array([0.01, 0.04, 0.03, 0.2])
    q = rs.bh(p)
    assert np.all(q >= p) and q.max() <= 1


def test_ora_finds_enriched_set():
    bg = [f"G{i}" for i in range(1000)]
    sets = {"HIT": [f"G{i}" for i in range(30)], "RANDOM": [f"G{i}" for i in range(500, 530)]}
    res = rs.ora([f"G{i}" for i in range(20)], bg, sets, 10, 500)
    assert res.iloc[0]["Term"] == "HIT" and res.iloc[0]["FDR"] < 0.001


def test_map_symbols_with_versions():
    ann = pd.DataFrame({"gene_name": ["TP53"]}, index=["ENSG00000141510"])
    out = rs.map_symbols(["ENSG00000141510.17", "ENSGX"], ann)
    assert out.tolist() == ["TP53", "ENSGX"]


def test_chromosome_harmonization(tmp_path, monkeypatch):
    gtf = tmp_path / "g.gtf"
    gtf.write_text('chr1\tx\texon\t1\t10\t.\t+\t.\tgene_id "A";\nchrM\tx\texon\t1\t10\t.\t+\t.\tgene_id "B";\n')
    monkeypatch.setattr(rs, "bam_contigs", lambda b: {"1", "MT"})
    fixed = rs.harmonize_gtf(str(gtf), ["x.bam"], tmp_path, lambda m: None)
    assert [l.split("\t")[0] for l in Path(fixed).read_text().splitlines()] == ["1", "MT"]


def test_nospace(tmp_path, monkeypatch):
    monkeypatch.setattr(rs, "APP_HOME", tmp_path / "home")
    d = tmp_path / "NUOVA ANALISI"
    d.mkdir()
    p = rs.nospace(d / "x.bam")
    assert " " not in p and Path(p).parent.resolve() == d.resolve()


@pytest.mark.slow
def test_deseq_end_to_end(tmp_path):
    counts, *_ = rs.read_counts(ROOT / "examples" / "example_counts.tsv")
    groups = {s: s.split("_")[0] for s in counts.columns if not s.startswith("TRT")}

    class _E:
        def __init__(self):
            self.emit = lambda *a: None

    class _T:
        log, progress = _E(), _E()

    cfg = dict(outdir=str(tmp_path), threads=2, alpha=0.05, lfc=1, min_count=10, min_samples=0, shrink=False,
               only_contrast=False, use_batch=False, annotation="", symbols=None)
    r = rs.task_deseq(_T(), counts, groups, {}, "WT", "MUT", cfg)
    up = r["res"].loc[[f"ENSG{i:011d}" for i in range(80)], "direction"]
    assert (up == "UP").mean() > 0.8
