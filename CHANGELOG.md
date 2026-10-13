# Changelog

## 1.0.0 — 2026-10-13

First public release (English interface).

- Four independent steps: alignment (fastp, HISAT2/STAR), counting (featureCounts), differential expression (PyDESeq2), pathway analysis (ORA, preranked GSEA).
- Self-installation without administrator rights: Python packages in `~/.local`, bioinformatics tools from Bioconda in `~/.local/share/rnaseqstudio/tools`, missing Qt system libraries from conda-forge.
- Per-sample roles (control, case, in model, excluded) assignable to single samples or whole groups; queue of multiple comparisons.
- Resumable runs, skipping of failed samples and corrupted BAMs, automatic `chr` prefix harmonisation between GTF and BAM, support for paths with spaces.
- Non-blocking warning and blind dispersion estimation when replicates are missing.
