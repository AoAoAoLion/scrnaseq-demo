#!/usr/bin/env python3
"""Single-cell RNA-seq analysis of a real 10x dataset (PBMC 3k).

Dataset: 10x Genomics PBMC 3k — 2,700 peripheral blood mononuclear cells from one donor,
released by 10x and distributed with scanpy. It is real data, publicly available, and small
enough that the whole analysis runs in a couple of minutes on a laptop.

Pipeline: QC -> normalisation -> HVGs -> PCA -> Leiden clustering -> UMAP -> marker genes ->
marker-based annotation -> HTML report.

Usage:
    python scripts/analyze.py --outdir results
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
from jinja2 import Template

sc.settings.verbosity = 1
plt.rcParams.update({"figure.dpi": 130, "font.size": 9})

# canonical PBMC markers: cluster identity is assigned only when the marker is both
# significantly enriched and expressed in a reasonable fraction of the cluster
MARKERS = {
    "T cells (CD4/CD8)": ["CD3D", "CD3E", "IL7R", "LTB"],
    "B cells": ["MS4A1", "CD79A", "CD79B"],
    "NK cells": ["NKG7", "GNLY", "GZMB"],
    "Monocytes (CD14)": ["LYZ", "S100A8", "S100A9", "CD14"],
    "Dendritic cells": ["FCER1A", "CST3"],
    "Platelets": ["PPBP", "PF4"],
    "Proliferating cells": ["MKI67", "TYMS", "PCNA", "TOP2A", "ZWINT"],
}


def log(msg: str) -> None:
    print(f"[scRNA] {msg}", flush=True)


def qc(adata, outdir: Path, min_genes: int, max_mito: float):
    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True)

    before = adata.n_obs
    fig, axes = plt.subplots(1, 3, figsize=(12, 3.4))
    axes[0].hist(adata.obs["n_genes_by_counts"], bins=60, color="#4C72B0")
    axes[0].axvline(min_genes, color="r", ls="--", lw=1)
    axes[0].set_xlabel("genes per cell")
    axes[1].hist(adata.obs["total_counts"], bins=60, color="#55A868")
    axes[1].set_xlabel("UMI counts per cell")
    axes[2].hist(adata.obs["pct_counts_mt"], bins=60, color="#C44E52")
    axes[2].axvline(max_mito, color="r", ls="--", lw=1)
    axes[2].set_xlabel("% mitochondrial counts")
    fig.suptitle(f"QC before filtering ({before:,} cells)")
    fig.tight_layout()
    fig.savefig(outdir / "figures" / "01_qc_before.png")
    plt.close(fig)

    sc.pp.filter_cells(adata, min_genes=min_genes)
    adata = adata[adata.obs["pct_counts_mt"] < max_mito].copy()
    sc.pp.filter_genes(adata, min_cells=3)
    log(f"QC: {before:,} -> {adata.n_obs:,} cells, {adata.n_vars:,} genes "
        f"(min_genes={min_genes}, max_mito={max_mito}%)")
    return adata, before


def cluster(adata, resolution: float, n_pcs: int, n_hvg: int):
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    adata.raw = adata
    sc.pp.highly_variable_genes(adata, n_top_genes=n_hvg, flavor="seurat")
    sc.pp.scale(adata, max_value=10)
    sc.tl.pca(adata, n_comps=n_pcs, svd_solver="arpack")
    sc.pp.neighbors(adata, n_neighbors=10, n_pcs=n_pcs)
    sc.tl.umap(adata)
    sc.tl.leiden(adata, resolution=resolution, key_added="cluster", flavor="igraph", n_iterations=2)
    log(f"clustering: {adata.obs['cluster'].nunique()} clusters at resolution {resolution}")
    return adata


def figures(adata, outdir: Path):
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.6))
    sc.pl.umap(adata, color="cluster", ax=axes[0], show=False, legend_loc="on data",
               legend_fontsize=7, frameon=False, title="Leiden clusters")
    sc.pl.umap(adata, color="total_counts", ax=axes[1], show=False, frameon=False,
               title="library size", cmap="viridis")
    fig.tight_layout()
    fig.savefig(outdir / "figures" / "02_umap.png")
    plt.close(fig)

    markers = ["CD3D", "MS4A1", "NKG7", "LYZ", "FCER1A", "PPBP"]
    present = [m for m in markers if m in adata.raw.var_names]
    fig, axes = plt.subplots(2, 3, figsize=(12, 6.5))
    for ax, gene in zip(axes.ravel(), present):
        sc.pl.umap(adata, color=gene, ax=ax, show=False, frameon=False, title=gene,
                   cmap="magma", vmin=0, vmax="p99")
    for ax in axes.ravel()[len(present):]:
        ax.axis("off")
    fig.tight_layout()
    fig.savefig(outdir / "figures" / "03_marker_umap.png")
    plt.close(fig)


def markers_and_annotation(adata, outdir: Path, min_frac: float = 0.25):
    sc.tl.rank_genes_groups(adata, "cluster", method="wilcoxon", use_raw=True)
    top = sc.get.rank_genes_groups_df(adata, group=None)
    top.to_csv(outdir / "tables" / "cluster_markers.csv", index=False)

    # score each cluster against the marker sets, weighted by detection rate
    rows = []
    for cluster in sorted(adata.obs["cluster"].unique(), key=int):
        sub = adata[adata.obs["cluster"] == cluster]
        best, best_score = "unassigned", 0.3   # below this, do not force a label
        scores = {}
        for label, genes in MARKERS.items():
            present = [g for g in genes if g in adata.raw.var_names]
            if not present:
                continue
            expr = np.asarray(sub.raw[:, present].X.toarray() if hasattr(
                sub.raw[:, present].X, "toarray") else sub.raw[:, present].X)
            frac = (expr > 0).mean(axis=0)
            score = float(np.mean(frac))
            scores[label] = round(score, 3)
            if score > best_score and (frac >= min_frac).sum() >= max(1, len(present) // 2):
                best, best_score = label, score
        rows.append({"cluster": cluster, "n_cells": int(sub.n_obs),
                     "assigned": best, "score": round(best_score, 3),
                     "scores": json.dumps(scores)})
    annot = pd.DataFrame(rows)
    annot.to_csv(outdir / "tables" / "cluster_annotation.csv", index=False)
    adata.obs["cell_type"] = [annot.set_index("cluster").loc[c, "assigned"]
                              for c in adata.obs["cluster"]]
    log("annotation: " + ", ".join(f"c{r.cluster}={r.assigned}" for r in annot.itertuples()))
    return annot


TEMPLATE = """<!doctype html><html lang="en"><head><meta charset="utf-8">
<title>scRNA-seq report — 10x PBMC 3k</title>
<style>
body{font-family:-apple-system,Segoe UI,Roboto,Helvetica,Arial,sans-serif;max-width:980px;margin:40px auto;padding:0 20px;color:#1b2430;line-height:1.6}
h1{border-bottom:3px solid #2ec4b6;padding-bottom:10px}
h2{margin-top:34px;border-bottom:1px solid #e3e8ef;padding-bottom:6px}
table{border-collapse:collapse;width:100%;font-size:14px;margin:10px 0}
th,td{border:1px solid #e3e8ef;padding:7px 9px;text-align:left}
th{background:#f6f8fb}
img{max-width:100%;border:1px solid #e3e8ef;border-radius:6px;margin:8px 0}
code{background:#f2f4f7;padding:2px 5px;border-radius:3px;font-size:13px}
.meta{color:#5b6b7c;font-size:14px}
</style></head><body>
<h1>Single-cell RNA-seq report — 10x PBMC 3k</h1>
<p class="meta">Dataset: 10x Genomics PBMC 3k (public) &middot; {{ n_cells }} cells after QC &middot;
{{ n_genes }} genes &middot; {{ n_clusters }} clusters &middot; scanpy {{ version }} &middot; {{ date }}</p>

<h2>1. Cells and QC</h2>
<table><tr><th>Step</th><th>Cells</th></tr>
<tr><td>Loaded</td><td>{{ n_before }}</td></tr>
<tr><td>After filter (min {{ min_genes }} genes, mito &lt; {{ max_mito }}%)</td><td>{{ n_cells }}</td></tr>
</table>
<img src="figures/01_qc_before.png" alt="QC distributions">
<p>Filtering thresholds are explicit and applied once. Cells outside them are removed before
normalisation, not after clustering.</p>

<h2>2. Clustering</h2>
<img src="figures/02_umap.png" alt="UMAP">
<p>{{ n_hvg }} highly variable genes, {{ n_pcs }} PCs, 10 neighbours, Leiden at resolution
{{ resolution }}. Clustering is on the expression graph; UMAP is used only for display.</p>

<h2>3. Cell types</h2>
<table><tr><th>Cluster</th><th>Cells</th><th>Assigned</th><th>Marker score</th></tr>
{% for r in annot %}<tr><td>{{ r.cluster }}</td><td>{{ r.n_cells }}</td><td>{{ r.assigned }}</td><td>{{ r.score }}</td></tr>{% endfor %}
</table>
<p>Labels come from canonical PBMC marker sets, scored by detection rate per cluster. A cluster
is left <code>unassigned</code> when no marker set passes both a score and a fraction-expressing
test, rather than being forced into the nearest label.</p>
<img src="figures/03_marker_umap.png" alt="markers">

<h2>4. Top markers per cluster</h2>
<table><tr><th>Cluster</th><th>Gene</th><th>log2FC</th><th>adj. p</th></tr>
{% for r in top_markers %}<tr><td>{{ r.group }}</td><td>{{ r.names }}</td><td>{{ '%.2f'|format(r.logfoldchanges) }}</td><td>{{ '%.1e'|format(r.pvals_adj) }}</td></tr>{% endfor %}
</table>
<p>Wilcoxon rank-sum on the raw counts, one test per cluster. Full table in
<code>results/tables/cluster_markers.csv</code>.</p>

<h2>5. Methods paragraph</h2>
<p>Cells with fewer than {{ min_genes }} detected genes or more than {{ max_mito }}% mitochondrial
counts were excluded. Counts were normalised to 10,000 per cell and log-transformed.
{{ n_hvg }} highly variable genes were selected, scaled, and reduced with PCA ({{ n_pcs }} components).
A k-nearest-neighbour graph (k = 10) was clustered with Leiden at resolution {{ resolution }}.
Cluster markers were identified with the Wilcoxon rank-sum test, and clusters were annotated
against canonical PBMC marker sets.</p>

<h2>6. Files</h2>
<ul>
<li><code>results/tables/cluster_markers.csv</code> — all markers, all clusters</li>
<li><code>results/tables/cluster_annotation.csv</code> — cluster sizes and assignment scores</li>
<li><code>results/adata_processed.h5ad</code> — processed object, re-openable with scanpy</li>
<li><code>results/figures/</code> — QC, UMAP, marker expression</li>
</ul>
</body></html>
"""


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--outdir", default="results")
    ap.add_argument("--min-genes", type=int, default=200)
    ap.add_argument("--max-mito", type=float, default=20.0)
    ap.add_argument("--n-hvg", type=int, default=2000)
    ap.add_argument("--n-pcs", type=int, default=50)
    ap.add_argument("--resolution", type=float, default=0.5)
    args = ap.parse_args()

    outdir = Path(args.outdir)
    (outdir / "figures").mkdir(parents=True, exist_ok=True)
    (outdir / "tables").mkdir(parents=True, exist_ok=True)

    log("loading 10x PBMC 3k (public dataset)")
    adata = sc.datasets.pbmc3k()
    adata.var_names_make_unique()

    adata, before = qc(adata, outdir, args.min_genes, args.max_mito)
    adata = cluster(adata, args.resolution, args.n_pcs, args.n_hvg)
    figures(adata, outdir)
    annot = markers_and_annotation(adata, outdir)

    top = pd.read_csv(outdir / "tables" / "cluster_markers.csv")
    top = top.sort_values(["group", "pvals_adj"]).groupby("group").head(3)

    ctx = dict(
        n_cells=int(adata.n_obs), n_genes=int(adata.n_vars), n_before=int(before),
        n_clusters=int(adata.obs["cluster"].nunique()), min_genes=args.min_genes,
        max_mito=args.max_mito, n_hvg=args.n_hvg, n_pcs=args.n_pcs,
        resolution=args.resolution, version=sc.__version__,
        date=pd.Timestamp.now().strftime("%Y-%m-%d %H:%M"),
        annot=annot.to_dict("records"), top_markers=top.to_dict("records"),
    )
    (outdir / "report.html").write_text(Template(TEMPLATE).render(**ctx))
    adata.write_h5ad(outdir / "adata_processed.h5ad")
    (outdir / "summary.json").write_text(json.dumps(
        {k: ctx[k] for k in ("n_cells", "n_genes", "n_clusters", "n_before")}, indent=2))
    log(f"done: {ctx['n_clusters']} clusters, report -> {outdir/'report.html'}")


if __name__ == "__main__":
    main()
