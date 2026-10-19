"""End-to-end test for the scRNA-seq demo.

Runs the real pipeline on the public PBMC 3k dataset with reduced parameters, so the test
exercises the actual code path (QC, clustering, markers, annotation, report) rather than a stub.
"""
from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

import pandas as pd
import pytest

REPO = Path(__file__).resolve().parents[1]
PY = sys.executable


@pytest.fixture(scope="module")
def run(tmp_path_factory):
    out = tmp_path_factory.mktemp("sc") / "results"
    proc = subprocess.run(
        [PY, str(REPO / "scripts" / "analyze.py"), "--outdir", str(out),
         "--n-hvg", "800", "--n-pcs", "20", "--resolution", "0.5"],
        capture_output=True, text=True, check=True,
    )
    return out, proc.stdout


def test_clusters_found(run):
    out, _ = run
    summary = json.loads((out / "summary.json").read_text())
    assert summary["n_cells"] > 2000
    assert summary["n_clusters"] >= 4


def test_marker_table_schema(run):
    out, _ = run
    markers = pd.read_csv(out / "tables" / "cluster_markers.csv")
    for col in ["group", "names", "logfoldchanges", "pvals_adj"]:
        assert col in markers.columns
    assert len(markers) > 500
    assert markers["pvals_adj"].notna().any()


def test_every_cluster_annotated(run):
    out, _ = run
    annot = pd.read_csv(out / "tables" / "cluster_annotation.csv")
    summary = json.loads((out / "summary.json").read_text())
    assert len(annot) == summary["n_clusters"]
    assert annot["assigned"].notna().all()
    assert annot["n_cells"].sum() == summary["n_cells"]
    # T cells dominate PBMCs; if this fails the annotation logic is broken, not the data
    assert "T cells" in annot.loc[annot["n_cells"].idxmax(), "assigned"]


def test_report_and_figures(run):
    out, _ = run
    html = (out / "report.html").read_text()
    assert "Single-cell RNA-seq report" in html
    assert "Methods paragraph" in html
    for fig in ["01_qc_before.png", "02_umap.png", "03_marker_umap.png"]:
        assert (out / "figures" / fig).exists(), fig
