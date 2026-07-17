# LUAD vs LUSC: CD8 and TAM Biology in Neoadjuvant Anti-PD-1 Treated NSCLC

## Research Question

How does CD8 and TAM biology (transcriptional activity, subtype heterogeneity, cell-cell communication) differ between lung adenocarcinoma (LUAD) and lung squamous cell carcinoma (LUSC) in patients treated with neoadjuvant anti-PD-1 therapy?

## Dataset

- **Source**: GSE243013 (Liu et al., 2025, *Cell*), single-cell RNA/TCR-seq of surgical tumor samples post-neoadjuvant chemo-immunotherapy in NSCLC.
- **Cohort**: 243 patients total (63 LUAD, 180 LUSC), verified from GEO metadata by unique `sampleID` and cross-checked against the `pathological_response` distribution (MPR n=45, non-MPR n=112, pCR n=85, unknown n=1).
- **Note on cohort size**: this is slightly higher than the 234 patients reported in the original publication's abstract. The source of this discrepancy was not resolved; it does not appear related to sample-level duplication or response-category mislabeling (both were checked directly against the metadata).
- **Metadata file**: `GSE243013_NSCLC_immune_scRNA_metadata.csv.gz`, cell-level resolution, key columns include `sampleID`, `cancer_type`, `pathological_response`, `major_cell_type`, `sub_cell_type`.

**EGFR status note:** Not available in public metadata as a dedicated field. 
No EGFR TKI recorded in the `targeted_therapy` field for either arm. This is consistent with the cohort's inclusion criteria as reported in the source publication (Liu et al., 2025): patients with EGFR mutation or ALK rearrangement were excluded by design, since these patients would be eligible for targeted therapies instead. The cohort can therefore be considered predominantly EGFR wild-type by construction, not merely by inference from comparable trials.

## Repository Structure

```
LUAD_LUSC_project/
├── Data/                                 # raw and processed input data
├── Objects/                              # saved Seurat/R objects
├── Scripts/
│   ├── Bloc0_DataPrep/
│   ├── Bloc1_QC/
│   ├── Bloc2_Norm_HVG_PCA/
│   ├── Bloc3_Harmony/
│   ├── Bloc4_Clustering_UMAP/
│   ├── Bloc5_Annotation_TME/
│   ├── Bloc6_CD8/
│   ├── Bloc7_TAM/
│   └── Bloc8_Communication/
└── Results/
    ├── Bloc0_DataPrep/{Figures,Tables}
    ├── Bloc1_QC/{Figures,Tables}
    ├── Bloc2_Norm_HVG_PCA/{Figures,Tables}
    ├── Bloc3_Harmony/{Figures,Tables}
    ├── Bloc4_Clustering_UMAP/{Figures,Tables}
    ├── Bloc5_Annotation_TME/{Figures,Tables}
    ├── Bloc6_CD8/{Figures,Tables}
    ├── Bloc7_TAM/{Figures,Tables}
    └── Bloc8_Communication/{Figures,Tables}
```

**Naming convention**: within each `Scripts/BlocX_*/` folder, scripts are numbered starting at `01`, following the execution order *within that block* (not a global pipeline count). Each script carries a histology suffix, `LUAD` or `LUSC`, run in parallel on both arms before the final comparison. Example: `01_LUAD_QC.R` and `01_LUSC_QC.R` in `Bloc0_QC`.

Results (figures and tables) mirror this structure exactly: same block folder, same script number, same histology suffix, so that a given script and its output can always be matched by name alone.

All paths are relative, managed via the `here` package (`here()` / `here("Data", ...)`), for portability across Windows, WSL, and Docker without modification.

## Pipeline

| Block | Content | Status |
|---|---|---|
| Bloc0_DataPrep | Matrix loading (BPCells), LUAD/LUSC separation, patient-level clinical metadata | Done |
| Bloc1_QC | Quality control plots and filtering (per-arm thresholds) | Done |
| Bloc2_Norm_HVG_PCA | Normalization, highly variable genes, PCA | Done |
| Bloc3_Harmony | Batch integration | Done |
| Bloc4_Clustering_UMAP | Clustering, UMAP visualization | Done |
| Bloc5_Annotation_TME | TME annotation (canonical markers as primary reference, SingleR fine labels as cross-check) | To do |
| Bloc6_CD8 | ProjecTILs CD8 subtyping, TCR clonality (scRepertoire), UCell signature scoring, CollecTRI regulon activity | To do |
| Bloc7_TAM | TAM subtyping, UCell signature scoring, CollecTRI regulon activity | To do |
| Bloc8_Communication | MultiNicheNet differential cell-cell communication (LUAD vs LUSC) | To do |

Each block's methodology, chosen parameters, and rationale will be documented in this README once the block is stabilized and just before its commit.

## Tools and Methods

| Tool | Role in the project | Reference | Block |
|---|---|---|---|
| BPCells | Disk-backed matrix loading, avoids RAM overload | Parikh et al., Satija Lab | Bloc0 |
| Harmony | Batch integration | Korsunsky et al., 2019, *Nature Methods* | Bloc3 |
| SingleR | TME annotation (fine-grained labels, cross-check against canonical markers) | Aran et al., 2019, *Nature Immunology* | Bloc5 |
| ProjecTILs | CD8 subtyping via reference atlas projection | Andreatta, Carmona et al., 2021, *Nature Communications* | Bloc6 |
| scRepertoire | TCR clonality analysis | Yang, Safina, Nguyen, Tuong, Borcherding, 2025, *PLOS Computational Biology* | Bloc6 |
| UCell | Gene signature scoring (CD8 and TAM) | Andreatta & Carmona, 2021, *Computational and Structural Biotechnology Journal* | Bloc6, Bloc7 |
| decoupleR + CollecTRI | Transcription factor activity inference (CD8 and TAM) | decoupleR: Badia-i-Mompel et al., 2022, *Bioinformatics Advances* — CollecTRI: Müller-Dott et al., 2023, *Nucleic Acids Research* | Bloc6, Bloc7 |
| MultiNicheNet | Differential cell-cell communication, multi-sample multi-condition design | Browaeys et al., 2023, *bioRxiv* | Bloc8 |

## Methodology

### Bloc0_DataPrep

- **Separation**: LUAD/LUSC split performed once (Script 00), producing reference ID lists (`sampleID`, `cellID`) reused by all downstream scripts. No re-filtering of `cancer_type` elsewhere in the pipeline.
- **Verification**: patient counts checked against previously validated values (63 LUAD, 180 LUSC) at each step; script halts on mismatch.
- **Matrix loading**: BPCells used for disk-backed loading, avoiding full in-memory load of the raw count matrix.
- **Clinical metadata**: one row per patient, deduplicated by `sampleID`. A known typo (`"unknowm"` instead of `"unknown"`) in `pathological_response` was kept as-is for traceability rather than silently corrected.
- **Response categories**: LUAD has 3 response categories (MPR, non-MPR, pCR). LUSC has 4, including one patient labeled `"unknowm"` (unclassified). This patient is retained in all data and score calculations, but excluded from MPR vs non-MPR statistical comparisons and main figures, same treatment as pCR (see ROADMAP.md).
- **Total cell count cross-check**: LUAD (336,685) + LUSC (918,064) = 1,254,749 cells, matching the full unfiltered metadata row count.

### Bloc1_QC 

- **Metrics**: nFeature_RNA, nCount_RNA, and percent.mt (mitochondrial content, `^MT-` pattern) computed per cell, after minimal Seurat object creation filters (`min.cells=3`, `min.features=200`).
- **Threshold selection**: thresholds based on the 99th percentile of each metric's distribution, computed separately for LUAD and LUSC (never copied across histologies), and uniform across all patients and pathological response groups within each histology. Differentiating thresholds by response group was explicitly avoided, since `pathological_response` is the primary comparison axis (MPR vs non-MPR); group-specific thresholds risk introducing an artificial bias into that comparison.
- **Thresholds applied**:
  - LUAD: nFeature_RNA [200–5000], nCount_RNA [1000–28000], percent.mt < 10%
  - LUSC: nFeature_RNA [200–4300], nCount_RNA [1000–26000], percent.mt < 10%
- **Verification**: cells removed per pathological response group checked post-filtering to confirm no group was disproportionately affected. LUAD: MPR -1.8%, non-MPR -2.7%, pCR -1.4%. LUSC: MPR -3.1%, non-MPR -3.0%, pCR -2.5%, unclassified ("unknowm") -2.5%. Comparable removal rates across groups in both histologies support the validity of uniform thresholds for the downstream MPR vs non-MPR comparison.
- **Object portability**: raw counts converted from BPCells (disk-backed) to standard in-memory sparse matrix (`dgCMatrix`) at the end of QC filtering, ensuring downstream objects remain portable across operating systems/environments, independent of the original BPCells storage path.

### Bloc2_Norm_HVG_PCA — Methodology

- **Normalization**: LogNormalize (scale factor 10,000), standard method.
- **Highly variable genes**: 3,000 HVGs selected (vst method), above the Seurat default of 2,000. Rationale verified against `major_cell_type` metadata: the dataset is CD45+ sorted (B cell, Myeloid cell, T/NK cell only; no epithelial population), so the increased HVG count is intended to better capture fine-grained heterogeneity within immune subpopulations (e.g. TAM subtypes, CD8 states), not epithelial diversity as originally assumed in the inherited script.
- **Cell cycle regression**: cell cycle phase scored (`CellCycleScoring`, Seurat's built-in S/G2M gene sets) and regressed out during scaling (`vars.to.regress = c("S.Score", "G2M.Score")`). Added after initial clustering (pre-regression) revealed a cluster driven entirely by proliferation markers (TYMS, RRM2, AURKB, UBE2C) mixing multiple cell lineages rather than reflecting a single biological identity.
- **Scaling**: applied to HVGs only, not the full gene set, to limit memory usage.
- **PCA**: 50 components computed; number of PCs to retain determined from per-PC variance explained, using the point where successive PCs contribute near-negligible marginal variance (<0.1%).
- **PC1 diagnostic**: correlation with nCount_RNA checked to rule out a purely technical origin for the dominant component (LUAD, post-regression: r=-0.50, PC1=30.4% variance); correlation with percent.mt negligible (r≈0).
- **PCs retained**: 20 for LUAD (post-regression), consistent with the pre-regression choice; margin beyond the observed elbow (~PC15-16) to preserve signal from rare subpopulations.

### Bloc3_Harmony

- **Batch correction**: Harmony (`group.by.vars = "sampleID"`), correcting for inter-patient technical variation, applied on the 20 PCs retained in Bloc2.
- **Verification**: patient-level dispersion (SD of per-patient mean coordinate) compared before/after Harmony as a proxy for integration quality. LUAD (post cell-cycle regression): 5.72→5.03 (-12%), an improvement over the pre-regression result (-7%), consistent with Harmony correcting more effectively on a PC1 signal less confounded by proliferation.

### Bloc4_Clustering_UMAP

- **Clustering**: `FindNeighbors`/`FindClusters` on the Harmony embedding, 20 dimensions, resolution 0.5. LUAD (post cell-cycle regression): 18 clusters, down from 19 pre-regression, consistent with the disappearance of the mixed-lineage proliferation cluster.
- **UMAP**: computed on the Harmony embedding, same 20 dimensions.
- **Integration quality check**: patient distribution visually inspected on the UMAP; well-mixed colors across all major cell islands, no patient-dominated isolated cluster identified.
- **Note**: LUSC has not yet been re-run with cell cycle regression at time of writing; LUSC section to be updated once completed.

### Bloc5_Annotation_TME

**LUAD**

- **Marker-based annotation**: `FindAllMarkers` (Wilcoxon test via `presto`), top10 markers per cluster reviewed manually against canonical immune markers (see MARKERS_CANONICAL_REFERENCE.md). 18 clusters (post cell-cycle regression) annotated into 17 homogeneous cell types; one cluster (cycling cells, mixed lineage) left unresolved by this method alone.
- **SingleR cross-check**: `MonacoImmuneData` reference (immune-only, appropriate for this CD45+ sorted dataset), `label.fine` as primary cross-check. Strong concordance with manual annotation; apparent disagreements traced to categories absent from the reference (TRM, exhausted CD4, ILC3, GC B cell, macrophage, mast cell — a peripheral-blood-derived reference lacks tissue-resident cell categories).
- **Mixed cluster resolution**: the unresolved cluster's CD8 cells (identified via SingleR label) were further characterized by a targeted `FindAllMarkers` comparison across their SingleR-assigned subgroups, without presupposing correspondence to existing CD8 categories (avoiding confirmation bias). This revealed three additional CD8 states not previously captured: PDCD1+ exhausted, TCF7+ naive/memory, and KLRC2+ NK-like terminal, all cycling.
- **Consensus**: manual annotation retained for homogeneous clusters; cell-level SingleR/targeted resolution applied only to the mixed cluster. 0 unannotated cells.
- **Statistical comparison**: global chi-squared test on the full composition table (MPR vs non-MPR; pCR retained through annotation but excluded from statistical comparisons, see Results for total CD8/TAM cell counts). Targeted proportion tests (`prop.test`) performed for B cell (sampling-artifact check) and for the CD8/TAM categories central to this study's hypothesis.

**LUSC**
- *Not yet completed. LUAD annotation completed first per project workflow (each histology processed fully through Bloc5-Bloc8 before starting the next)*



## Preliminary Observations 

### Bloc5_Annotation_TME

**LUAD**
- 18 clusters identified post cell-cycle regression, annotated into 20 immune populations (17 homogeneous + 3 additional cycling CD8 substates resolved from one mixed cluster).
- TME composition differs significantly between MPR and non-MPR (chi-squared, p < 2.2e-16).
- B cell proportion markedly higher in MPR (41.2% vs 18.1%), consistent with a distinct GC B cell cluster suggesting possible tertiary lymphoid structures.
- Macrophage/TAM (C1Q+/lipid-associated) proportion nearly 3x higher in non-MPR (12.2% vs 4.6%), the largest effect size among tested categories.
- CD8 T (LAG3+ exhausted/TRM) also higher in non-MPR (12.4% vs 8.8%); CD8 T (GZMK+ effector/memory) shows a statistically significant but modest difference (12.4% vs 11.4%).
- Total CD8 T cells across all states: 80,615 cells, spanning 5 distinct populations.



## Discussion

*To be written once results are available. Condensed prose, distinct from block-level preliminary observations.*