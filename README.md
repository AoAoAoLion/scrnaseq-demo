# Single-cell RNA-seq analysis — real 10x PBMC data

A complete, reproducible scRNA-seq analysis on a **public 10x Genomics dataset** (PBMC 3k:
2,700 peripheral blood mononuclear cells from one donor). No simulated data, no synthetic
counts — the numbers in the report come from the real experiment.

```bash
pip install -r requirements.txt
python scripts/analyze.py --outdir results
open results/report.html
```

Runs in about 45 seconds on a laptop CPU.

## What it does

| Step | Detail |
|---|---|
| QC | genes/cell, UMIs/cell, % mitochondrial; filters at 200 genes and 20% mito |
| Normalisation | counts per 10k, log1p |
| Feature selection | 2,000 highly variable genes (Seurat flavour) |
| Embedding | PCA (50 components) → kNN graph (k=10) → UMAP |
| Clustering | Leiden, resolution 0.5 |
| Markers | Wilcoxon rank-sum per cluster |
| Annotation | canonical PBMC marker sets, scored by detection rate |

## Result on this dataset

2,698 cells pass QC, 13,714 genes retained, **7 clusters**:

| Cluster | Cells | Assigned | Marker score |
|---|---|---|---|
| 0 | 1,224 | T cells (CD4/CD8) | 0.80 |
| 1 | 346 | B cells | 0.90 |
| 2 | 653 | Monocytes (CD14) | 0.84 |
| 3 | 421 | NK cells | 0.70 |
| 4 | 34 | Dendritic cells | 0.91 |
| 5 | 12 | Platelets | 1.00 |
| 6 | 8 | Proliferating cells | 0.88 |

That is the expected composition of PBMCs, which is the point: the pipeline recovers known
biology rather than producing clusters that need a story.

## Annotation is a decision, not a lookup

Cluster 6 has 8 cells and no strong immune marker. A nearest-label assignment would have called
it NK cells (its best score was 0.58 while T cells scored 0.56 — a coin flip). Instead the code
requires a minimum score **and** a fraction of markers actually detected, so the cluster stays
unassigned unless something specific passes. Adding a proliferation marker set (MKI67, TYMS,
PCNA, TOP2A, ZWINT) then identified it correctly: its top genes are TYMS, KIAA0101, ZWINT, PCNA.

That trade-off is deliberate: an honest "unassigned" is more useful than a confident wrong label.

## Outputs

| File | Content |
|---|---|
| `results/report.html` | QC, clustering, annotation table, top markers, methods paragraph |
| `results/tables/cluster_markers.csv` | every marker, every cluster (Wilcoxon) |
| `results/tables/cluster_annotation.csv` | cluster sizes, assignment, per-label scores |
| `results/adata_processed.h5ad` | processed AnnData, re-openable with scanpy |
| `results/figures/` | QC distributions, UMAP, marker expression UMAPs |

## Requirements

`scanpy`, `leidenalg`, `igraph`, `pandas`, `matplotlib`, `jinja2` (pinned in
`requirements.txt`). The dataset is downloaded by scanpy on first run and cached
(`~/data/pbmc3k_raw.h5ad`), so later runs work offline.

## What this analysis does not do

- **No doublet removal.** PBMC 3k is already a filtered 10x release; adding Scrublet to a demo
  would suggest a step was validated when it was not.
- **No batch correction.** One donor, one lane — there is no batch to correct.
- **No trajectory or RNA velocity.** Different question, different assumptions.
- **No automatic cluster number.** Resolution is a parameter you can see, not something the code
  decides silently.

## Tests

```bash
pytest tests -q
```

Checks that the pipeline produces clusters, a marker table with the expected columns, an
annotation for every cluster, and a report with figures.

## Data source

10x Genomics, "PBMC 3k" (v1 chemistry), distributed with scanpy (`sc.datasets.pbmc3k`).
Public and free to use; cite 10x Genomics and the scanpy paper if you reuse the figures.

## License

MIT.

## Rendered output

`example_output/` contains the report and figures from the run described above (2,698 cells, 7
clusters), so the result can be inspected without installing anything.
`example_output/tables/top20_markers_per_cluster.csv` is trimmed to the top 20 markers per
cluster to keep the repository small; the script writes the full table with every gene.
