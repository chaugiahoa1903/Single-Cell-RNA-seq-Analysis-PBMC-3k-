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
   (at 0.2, the two CD4 T clusters merge into one, NK and CD8 T cells merge, and the two monocyte groups are not split);
   at 0.8 and 1.2 the T/NK region breaks into many small groups that I did not test further.
5. **Annotation.** Marker genes per cluster (Wilcoxon rank-sum) checked on a dotplot:

| Cell type | Main markers seen | % of cells |
|---|---|---|
| CD4 T (2 clusters, see item 7) | IL7R, LTB, CD3D | 43.6 |
| CD14+ monocyte | LYZ, S100A8/9, CD14 | 16.2 |
| B | MS4A1, CD79A, CD74 | 12.6 |
| CD8 T | CCL5, NKG7, GZMA, CD8A/B (weak) | 12.1 |
| FCGR3A+ monocyte | FCGR3A, MS4A7, LST1, FCER1G | 7.7 |
| NK | GNLY, NKG7, GZMB, PRF1 | 6.0 |
| Dendritic | FCER1A, CST3, HLA-DPA1 | 1.2 |
| Megakaryocyte | PPBP, PF4, NRGN | 0.6 |

6. **Differential expression: FCGR3A+ vs CD14+ monocytes.** Top genes: FCGR3A (log fold change 6.1, adjusted p = 4e-64),
   IFITM2, LST1, FCER1G, MS4A7. This matches the known profile of CD16+ (non-classical) monocytes.
7. **CD4 T sub-clusters (731 and 401 cells).**
   - *Quality check:* the smaller cluster (1) has fewer genes (median 694 vs 868) and counts (2,019 vs 2,508); mitochondrial
     content is similar (1.9% vs 1.7%).
   - *Differential expression:* cluster 0 is higher in S100A4 and ANXA1 (log fold change about -1.2 for cluster 1 vs 0),
     S100A11, SRGN, VIM and LTB; cluster 1 is higher in MALAT1 and ribosomal genes, with small fold changes (0.15-0.4).
   - *Memory vs naive markers:* the memory score is higher in cluster 0 (mean 0.74 vs 0.32). Naive markers (CCR7, SELL,
     TCF7, LEF1) are similar in both clusters (naive score -0.59 vs -0.43), so I call cluster 0 **memory-like** and
     did **not** label cluster 1 as naive.
   - *Depth-matched check:* among cells with 1,500-3,000 counts (573 and 298 cells), the memory score difference remains
     (median 0.80 vs 0.29), so the split is not only a sequencing-depth effect.
   - *Resolution:* the two clusters are one cluster at resolution 0.2 (731 + 397 cells) and separate at 0.5.
8. **Stability.** Re-running neighbors + Leiden with 5 other seeds always gave 9 clusters, with adjusted Rand index
   **0.76-0.93 (mean 0.81)** against the reference. The number of clusters is stable; boundaries between neighboring
   clusters move somewhat. For the CD4 T cells (1,132), the same two-way split was recovered
   for **79.5-95.9% of cells across the 5 seeds (mean 84.7%)**: seed 1 (ARI 0.93, the most stable) 95.9%, seeds 2-5 79.5-82.8%.
   So the CD4 split is reproducible in direction but its boundary moves for roughly 1 in 6 cells.

9. **Cross-check in R (Seurat).** Same data, thresholds (200-2,500 genes, < 5% mito), 10 PCs and resolution 0.5
   (`seurat_pbmc3k.R`). Seurat found **9 clusters / the same 8 cell types** (2,638 cells; it did not remove doublets).
   On the 2,595 cells shared with the Scanpy analysis:
   - Cell-type labels agree for **96.3%** of cells (adjusted Rand index **0.94**, CD4 naive + memory merged into CD4 T).
     Most disagreements are FCGR3A+ vs CD14+ monocytes (45 cells) and CD8 T vs CD4 T (31 cells).
   - Cluster-level ARI (9 vs 9 clusters) is lower, **0.66**, mainly because of how the CD4 T cells are split.
   - **CD4 T split:** Seurat separates *naive* (CCR7, LEF1, TCF7 high; 673 cells) from *memory* (IL7R, CD40LG, AQP3; 455 cells).
     Of Scanpy's lower-memory cluster (401 cells), 373 (93%) are Seurat naive; of Seurat memory cells, 428 of 455 (94%) are in
     Scanpy's memory-like cluster. But Scanpy's memory-like cluster also holds 300 Seurat-naive cells, so the two tools agree
     on the direction of the split for about 71% of CD4 T cells, not on its exact boundary.

## Limitations
- PBMC 3k is a small, well-known benchmark; the cell types found are expected, not novel.
- Annotation uses a small set of canonical markers plus manual review. FCGR3A is also expressed by NK cells, so it
  is not specific to monocytes on its own, and CD8A/CD8B signal is weak in this sparse dataset.
- The CD4 T split is memory-like vs lower-memory, based on a small marker set; naive markers were not clearly enriched,
  and cluster 1 has lower depth, so part of the split may still be technical. The depth check only controls total counts,
  and its boundary is not fixed: across 5 seeds, 79.5-95.9% of CD4 T cells keep the same split.
- The CD4 T split is only partly reproduced across tools (about 71% of CD4 T cells); treat it as a tentative naive-like vs
  memory-like gradient, not two sharply separated populations (they form one connected region on the UMAP).
- The Seurat run did not remove doublets, so cell counts differ slightly between the two analyses.
- Differential-expression p-values are optimistic because the groups were defined from the same data.
- Scrublet was run with default settings and its expected doublet rate was not tuned.
