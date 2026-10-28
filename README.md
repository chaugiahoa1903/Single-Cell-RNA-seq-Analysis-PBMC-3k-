# PBMC 3k single-cell RNA-seq analysis

Scanpy workflow on the public 10x Genomics PBMC 3k dataset (2,700 peripheral blood mononuclear cells,
32,738 genes): QC, doublet removal, clustering at several resolutions, marker-gene annotation,
differential expression, and a stability check.

## How to run
1. `pip install scanpy igraph scikit-learn scikit-image`
2. Open `pbmc3k_analysis.ipynb` (Jupyter or Google Colab) and run all cells.
   The dataset is downloaded automatically by Scanpy. Random seed is fixed (`SEED = 0`).
   Tested with Python 3.11 and Scanpy 1.11.5.

## What I did and what I found

1. **QC.** Kept cells with 200-2,500 detected genes and < 5% mitochondrial reads; kept genes seen in >= 3 cells.
   Most cells have roughly 500-1,500 genes and low mitochondrial content, with a thin tail of cells above
   2,500 genes (possible doublets) and above 5% mitochondrial reads (possibly dying cells), which the filters remove.
   **2,700 -> 2,638 cells, 13,714 genes.**
2. **Doublets.** Scrublet (built into Scanpy), automatic threshold 0.24: **43 predicted doublets removed -> 2,595 cells.**
   Scrublet's estimated doublet rate (4.3%) was close to its default expected rate (5%), but only 1.6% of cells were flagged.
3. **Preprocessing.** Library-size normalization (10,000 counts/cell), log1p, **1,848 highly variable genes**, scaling, PCA.
   Used **10 PCs**: the variance curve drops sharply by about PC 6-8 and is flat after about PC 10.
4. **Clustering.** Leiden at four resolutions gave **6 / 9 / 11 / 15 clusters** (0.2 / 0.5 / 0.8 / 1.2).
   I chose **0.5 (9 clusters)** because it is the lowest resolution that separates all the expected PBMC populations
   (at 0.2, NK and CD8 T cells merge and the monocyte groups are not split); at 0.8 and 1.2 the T/NK region breaks into many
   small groups that I did not test further.
5. **Annotation.** Marker genes per cluster (Wilcoxon rank-sum) checked on a dotplot:

| Cell type | Main markers seen | % of cells |
|---|---|---|
| CD4 T (2 clusters) | IL7R, LTB, CD3D (CCR7 higher in the second cluster) | 43.6 |
| CD14+ monocyte | LYZ, S100A8/9, CD14 | 16.2 |
| B | MS4A1, CD79A, CD74 | 12.6 |
| CD8 T | CCL5, NKG7, GZMA, CD8A/B (weak) | 12.1 |
| FCGR3A+ monocyte | FCGR3A, MS4A7, LST1, FCER1G | 7.7 |
| NK | GNLY, NKG7, GZMB, PRF1 | 6.0 |
| Dendritic | FCER1A, CST3, HLA-DPA1 | 1.2 |
| Megakaryocyte | PPBP, PF4, NRGN | 0.6 |

6. **Differential expression: FCGR3A+ vs CD14+ monocytes.** Top genes: FCGR3A (log fold change 6.1, adjusted p = 4e-64),
   IFITM2, LST1, FCER1G, MS4A7. This matches the known profile of CD16+ (non-classical) monocytes.
7. **Stability.** Re-running neighbors + Leiden with 5 other seeds always gave 9 clusters, with adjusted Rand index
   **0.76-0.93 (mean 0.81)** against the reference. The number of clusters is stable; boundaries between neighboring
   clusters move somewhat.

## Limitations
- PBMC 3k is a small, well-known benchmark; the cell types found are expected, not novel.
- Annotation uses a small set of canonical markers plus manual review. FCGR3A is also expressed by NK cells, so it
  is not specific to monocytes on its own, and CD8A/CD8B signal is weak in this sparse dataset.
- The two CD4 T clusters were not characterized further (for example naive vs memory), so I did not label them separately.
- Differential-expression p-values are optimistic because the groups were defined from the same data.
- Scrublet was run with default settings and its expected doublet rate was not tuned.
