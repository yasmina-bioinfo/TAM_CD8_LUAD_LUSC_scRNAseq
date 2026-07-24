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

---

## Pipeline

| Block | Content | Status |
|---|---|---|
| Bloc0_DataPrep | Matrix loading (BPCells), LUAD/LUSC separation, patient-level clinical metadata | Done |
| Bloc1_QC | Quality control plots and filtering (per-arm thresholds) | Done |
| Bloc2_Norm_HVG_PCA | Normalization, highly variable genes, PCA | Done |
| Bloc3_Harmony | Batch integration | Done |
| Bloc4_Clustering_UMAP | Clustering, UMAP visualization | Done |
| Bloc5_Annotation_TME | TME annotation (CellTypist as primary cell-level method, canonical markers as initial exploratory pass) | Done |
| Bloc6_CD8 | ProjecTILs CD8 subtyping, TCR clonality (scRepertoire), UCell signature scoring, CollecTRI regulon activity | To do |
| Bloc7_TAM | TAM subtyping, UCell signature scoring, CollecTRI regulon activity | To do |
| Bloc8_Communication | MultiNicheNet differential cell-cell communication (LUAD vs LUSC) | To do |

Each block's methodology, chosen parameters, and rationale will be documented in this README once the block is stabilized and just before its commit.

---

## Tools and Methods

| Tool | Role in the project | Reference | Block |
|---|---|---|---|
| BPCells | Disk-backed matrix loading, avoids RAM overload | Parikh et al., Satija Lab | Bloc0 |
| Harmony | Batch integration | Korsunsky et al., 2019, *Nature Methods* | Bloc3 |
| ProjecTILs | CD8 subtyping via reference atlas projection | Andreatta, Carmona et al., 2021, *Nature Communications* | Bloc6 |
| scRepertoire | TCR clonality analysis | Yang, Safina, Nguyen, Tuong, Borcherding, 2025, *PLOS Computational Biology* | Bloc6 |
| UCell | Gene signature scoring (CD8 and TAM) | Andreatta & Carmona, 2021, *Computational and Structural Biotechnology Journal* | Bloc6, Bloc7 |
| decoupleR + CollecTRI | Transcription factor activity inference (CD8 and TAM) | decoupleR: Badia-i-Mompel et al., 2022, *Bioinformatics Advances* - CollecTRI: Müller-Dott et al., 2023, *Nucleic Acids Research* | Bloc6, Bloc7 |
| MultiNicheNet | Differential cell-cell communication, multi-sample multi-condition design | Browaeys et al., 2023, *bioRxiv* | Bloc8 |
| SingleR | *Initially used for TME annotation cross-check; replaced by CellTypist after insufficient robustness for closely related T-cell states (see Methodology)* | Aran et al., 2019, *Nature Immunology* | Bloc5 |
| CellTypist | TME annotation (cell-level classification, majority-voting consensus), replaces SingleR | Domínguez Conde et al., 2022, *Science* | Bloc5 |

---

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
  - LUAD: nFeature_RNA [200-5000], nCount_RNA [1000-28000], percent.mt < 10%
  - LUSC: nFeature_RNA [200-4300], nCount_RNA [1000-26000], percent.mt < 10%
- **Verification**: cells removed per pathological response group checked post-filtering to confirm no group was disproportionately affected. LUAD: MPR -1.8%, non-MPR -2.7%, pCR -1.4%. LUSC: MPR -3.1%, non-MPR -3.0%, pCR -2.5%, unclassified ("unknowm") -2.5%. Comparable removal rates across groups in both histologies support the validity of uniform thresholds for the downstream MPR vs non-MPR comparison.
- **Object portability**: raw counts converted from BPCells (disk-backed) to standard in-memory sparse matrix (`dgCMatrix`) at the end of QC filtering, ensuring downstream objects remain portable across operating systems/environments, independent of the original BPCells storage path.
- **Response group filtering**: following QC filtering, cells are restricted to `pathological_response %in% c("MPR", "non-MPR")`. pCR patients (and the single "unknowm"-labeled LUSC patient) are excluded from this point onward in both histologies, since neither is used in any downstream scoring or comparison beyond a descriptive pCR/total ratio per histology (reported separately from clinical metadata, not from the filtered Seurat object). Cell counts after filtering: LUAD 268,615 cells (MPR 50,914 / non-MPR 217,701, from 328,854 before filtering); LUSC 478,223 cells (MPR 183,956 / non-MPR 294,267, from 892,413 before filtering).

### Bloc2_Norm_HVG_PCA

- **Normalization**: LogNormalize (scale factor 10,000), standard method.
- **Highly variable genes**: 3,000 HVGs for both LUAD and LUSC (vst method), above the Seurat default of 2,000. Verified against `major_cell_type` metadata (CD45+ sorted dataset, no epithelial compartment), intended to better capture fine-grained heterogeneity within immune subpopulations. LUSC initially used a reduced HVG count (2,000) as a temporary computational constraint on the larger pre-filtering dataset; unified to 3,000 (matching LUAD) once the response-group filtering (see Bloc1_QC) reduced the LUSC cell count enough to remove this constraint.
- **Cell cycle regression**: not applied to either histology. Initially applied to LUAD only, then reassessed: regression reduced but did not fully resolve a mixed-lineage proliferation cluster (a residual mixed cluster still required post-hoc resolution via SingleR - see Bloc5), so it did not provide a complete solution even where it succeeded technically. On LUSC, regression additionally caused repeated memory allocation failures on the larger dataset. Since post-hoc SingleR resolution alone proved sufficient on LUAD, regression was dropped entirely for consistency across both histologies.
- **Scaling**: applied to HVGs only, not the full gene set, to limit memory usage. The dense scaled matrix (`scale.data`) is discarded from the object immediately after PCA, since it is not required by any downstream step (Harmony, clustering, UMAP, `FindAllMarkers`, SingleR all operate on the PCA/Harmony embedding or on the normalized sparse matrix).
- **PCA**: 50 components computed for both arms; number of PCs to retain determined from per-PC variance explained, using the point where successive PCs contribute a near-negligible marginal variance (<0.1%).
- **PC1 diagnostic**: correlation with nCount_RNA checked to rule out a purely technical origin for the dominant component (LUAD: r=-0.50, 30.8% variance; LUSC: r=0.48, 28.9% variance - both moderate, not indicative of a purely technical artifact); correlation with percent.mt negligible in both arms (LUAD r≈-0.02; LUSC r≈-0.06).
- **PCs retained**: 20 for both LUAD and LUSC, slightly beyond the strict elbow (~PC12-13) to reduce the risk of losing signal from rare subpopulations of interest (TAM, CD8 substates). Cumulative 20-PC variance: 82.4% (LUAD), 82.4% (LUSC).

### Bloc3_Harmony

- **Batch correction**: Harmony (`group.by.vars = "sampleID"`), correcting for inter-patient technical variation, applied on the 20 PCs retained in Bloc2.
- **Verification**: patient-level dispersion (SD of per-patient mean coordinate) compared before/after Harmony as a proxy for integration quality. LUAD: 5.93→5.40 (-8.9%). LUSC: 5.53→5.03 (-9.0%). Comparable reductions across both histologies support the robustness of the integration step.

### Bloc4_Clustering_UMAP

- **Clustering**: `FindNeighbors`/`FindClusters` on the Harmony embedding (not raw PCA), 20 dimensions, resolution 0.5 as a conservative starting point for global TME clustering (17 clusters LUAD; 21 clusters LUSC).
- **UMAP**: computed on the Harmony embedding, same 20 dimensions.
- **Integration quality check**: patient distribution visually inspected on the UMAP; well-mixed colors across all major cell islands in both LUAD and LUSC.

### Bloc5_Annotation_TME

**LUAD**

- **Marker-based annotation**: `FindAllMarkers` (Wilcoxon test via `presto`), top10 markers per cluster reviewed manually against canonical immune markers. 17 clusters (filtered MPR/non-MPR object, no cell-cycle regression) reviewed as an initial exploratory pass; not merged into the final annotation (see below).
- **Cell-level annotation (CellTypist)**: `Immune_All_Low` model (Domínguez Conde et al., 2022, *Science*; multi-tissue atlas including tumor-infiltrating lymphocytes), run with `majority_voting = True`. Applied systematically to all 17 clusters; a cluster was confirmed directly when ≥95% of its cells shared a CellTypist label within a single biological family. This directly confirmed 12/17 clusters. Direct comparison against an initial SingleR (MonacoImmuneData) cross-check on the same 17 clusters showed CellTypist achieving broader and cleaner direct confirmation (12/17 vs. 11/17 clusters confirmed directly). Mixed clusters requiring cell-level resolution also showed markedly fewer sub-groups with CellTypist: after filtering out statistically unreliable labels, 8-12 sub-groups per cluster (CellTypist) versus up to 14-15 (SingleR), out of 8-12 versus up to 22-25 raw distinct labels respectively. This is consistent with SingleR/MonacoImmuneData's peripheral-blood origin being less suited to tissue-infiltrating and closely related T-cell states. CellTypist was therefore adopted as the sole cell-level annotation method.
- **Mixed cluster resolution**: for the 5 remaining clusters (1, 2, 3, 6, 11), cells were grouped by their individual `CellTypist_majority` label, then a targeted `FindAllMarkers` was run per subgroup. Each subgroup's marker profile was checked against the CellTypist authors' own curated marker sets (Mullan et al., 2023, *Front Immunol* for T-cell subtypes; Domínguez Conde et al., 2022 for NK/myeloid/B) — a subgroup was confirmed only if its official marker(s) were recovered among its own top differentially expressed genes; unconfirmed subgroups were retained as explicitly traceable "Unresolved" labels rather than assigned an unverified identity. Signals dominated by non-immune markers (lung epithelial contamination: *SCGB1A1*, *SFTPB*, *SFTPC*, *SCGB3A1*) were excluded as non-immune contamination — a pattern independently reported by the CellTypist authors as an ambient-RNA/doublet artifact in their own macrophage cluster analysis.
- **Final annotation** (`annotation_final`): built entirely from CellTypist. Each retained label is the model's own terminology, not translated into project-specific nomenclature, to preserve direct traceability to the source classifier. 268,615 cells, 0 unannotated. 86.5% of cells received a marker-confirmed identity; 13.0% remained "Unresolved" and 0.6% were excluded as non-immune contamination. Technical lineage flags (`lineage_CD8`, `lineage_TAM`, `lineage_TAM_related`) were built from confirmed labels only, for downstream Bloc6/Bloc7 extraction.
- **Display simplification**: for figure readability only, categories below 500 cells are grouped into a single "Minor/unresolved categories" bin, and non-immune contamination into a separate bin; this applies to the summary figure (`annotation_display`) only and does not affect `annotation_final`, which retains every individually confirmed or unconfirmed label and is used for all statistical tests and downstream analyses.
- **Statistical comparison**: global chi-squared test on the full composition table (MPR vs non-MPR; pCR and the single "unknowm"-labeled patient are not present in the object at this stage). Targeted proportion tests (`prop.test`) performed for B cells (sampling-artifact check) and for the CD8/TAM categories central to this study's hypothesis, computed both per individual subtype and pooled across all marker-confirmed CD8/TAM states.

**LUSC**

- **Marker-based annotation**: `FindAllMarkers` (Wilcoxon test via `presto`), top10 markers per cluster reviewed manually against canonical immune markers. 21 clusters (filtered MPR/non-MPR object, no cell-cycle regression) reviewed as an initial exploratory pass; not merged into the final annotation.
- **Patient composition check**: during marker-based cluster review, one myeloid cluster (391 cells) showed an unexpected *XIST* signal in its top differentially expressed genes — a sex-linked marker with no lineage-specific relevance, prompting a check of patient composition. This cluster was found to originate entirely from a single patient (100%); the patient's overall cell count and QC metrics were within normal range across other clusters, and the pattern was attributed to a residual Harmony integration artifact for this specific cell subset rather than a quality issue. The patient was retained in the analysis without modification; CellTypist confirmed a biologically plausible identity (Classical monocytes) for this cluster.
- **Cell-level annotation (CellTypist)**: `Immune_All_Low` model (Domínguez Conde et al., 2022, *Science*), run with `majority_voting = True`, applied identically to LUSC as to LUAD. Applied systematically to all 21 clusters using the same ≥95% single-family concordance threshold. This directly confirmed 13/21 clusters; the remaining 8 (clusters 1, 2, 3, 4, 6, 12, 14, 20) mixed multiple cell states and required resolution at the single-cell level.
- **Mixed cluster resolution**: for each of the 8 mixed clusters, cells were grouped by their individual `CellTypist_majority` label, then a targeted `FindAllMarkers` was run per subgroup. Each subgroup's marker profile was checked against the CellTypist authors' own curated marker sets (Mullan et al., 2023; Domínguez Conde et al., 2022) — a subgroup was confirmed only if its official marker(s) were recovered among its own top differentially expressed genes; unconfirmed subgroups were retained as explicitly traceable "Unresolved" labels. A recurrent lung epithelial contamination signal (*SCGB1A1*, *SFTPC*, *SFTPA1*, *SFTPB*) was identified within otherwise immune-labeled subgroups in several clusters and excluded, consistent with the pattern observed in LUAD and with the ambient-RNA/doublet artifact independently reported by the CellTypist authors.
- **Final annotation** (`annotation_final`): built entirely from CellTypist, using the model's own terminology without translation to project-specific nomenclature. 478,223 cells, 0 unannotated. One resolved cluster (391 cells, confirmed as Classical monocytes) originated entirely from a single patient and is flagged as a Harmony integration artifact rather than a distinct biological subpopulation, though its cell-level identity remains biologically plausible. Technical lineage flags (`lineage_CD8`, `lineage_TAM`, `lineage_TAM_related`) were built from marker-confirmed labels only, for downstream Bloc6/Bloc7 extraction; LUSC additionally distinguishes "Macrophages" from "Alveolar macrophages" as two separately confirmed labels, both included in `lineage_TAM`.
- **Display simplification and statistical comparison**: identical approach to LUAD (see above) — display-only grouping of small/unresolved categories, global chi-squared test on the full composition table, and targeted `prop.test` for B cells, CD8, and TAM categories (individual subtypes and pooled marker-confirmed states).

### Bloc6_CD8

**LUAD**

- **CD8 purification strategy**: rather than subsetting on the CellTypist-derived `lineage_CD8` flag, CD8 T cells were re-isolated independently on the full TME object using scGate (`TME_HiRes` model, human), following the same isolation-then-projection strategy used in the official construction of the human CD8 TIL reference atlas (Andreatta et al., 2021, *Nat Commun*; atlas construction pipeline: scGate → UCell cycling/ISG+ removal → STACAS). This avoids circularity with the project's own CellTypist-based annotation and keeps downstream steps methodologically consistent with the reference space used for later ProjecTILs projection.
- **Performance note**: scGate recalculates a dimensionality reduction at each hierarchical gating level by default; supplying the pre-computed Harmony embedding (`reduction = "harmony"`) instead avoids this recomputation and was necessary for tractable runtime (268,615 cells, ~1h17 with this optimization).
- **scGate_multi distribution**: Bcell 60,743; CD4T 72,008; CD8T 58,203; Macrophage 18,220; Monocyte 8,876; NK 7,642; panDC 7,709; PlasmaCell 7,237; Mast 3,982; Neutrophils 1,486; Multi (ambiguous) 3,591; unassigned 18,918. Endothelial, Epithelial, and Fibroblast levels returned 0 cells each, consistent with the CD45+-sorted nature of the dataset.
- **Cross-check against CellTypist**: of the 57,782 cells confirmed as `lineage_CD8` by CellTypist, 41,752 (72.3%) were independently confirmed as pure CD8T by scGate. The 27.7% discordance reflects the different underlying approaches (generic marker-based gating vs. reference-based classification) rather than a methodological error; biological interpretation of this discordance is deferred until after ProjecTILs projection, which uses a CD8-specific reference atlas. `scGate_multi == "CD8T"` (not `lineage_CD8`, and not an intersection/union of the two) was retained as the base population carried forward to the next step, to preserve consistency with the reference atlas construction protocol; `lineage_CD8` is retained as an annotation/cross-validation column only.
- **Cycling/IFN-high cell removal (UCell)**: applied to the 58,203 scGate-confirmed CD8T cells, using the identical signatures and thresholds from the official CD8 human TIL atlas construction pipeline (`cellCycle.G1S`, `cellCycle.G2M` from SignatuR, Tirosh et al., 2016, *Science*; a manually-defined 4-gene IFN signature — *ISG15, IFI6, IFI44L, MX1* — reported by the atlas authors to perform better with UCell than the full SignatuR IFN signature). Thresholds: `cellCycle.G1S < 0.1`, `cellCycle.G2M < 0.1`, `IFN < 0.25`. 8,616 cells (14.8%) were removed; the IFN criterion accounted for the large majority of flagged cells (6,974, ~81% of cumulative flags) compared to the cell-cycle criteria combined (1,543 + 1,392). 49,587 cells retained for the next step (ProjecTILs projection).
- **ProjecTILs functional-state projection**: the 49,587 scGate/UCell-purified CD8 cells were projected onto the official human CD8 TIL reference atlas (Andreatta et al., 2021, *Nat Commun*; atlas construction: `carmonalab/CD8_human_TIL_atlas_construction`), using `Run.ProjecTILs()` with `filter.cells = FALSE` (query already pre-purified, per official guidance for pre-filtered data) and `split.by = "sampleID"` (per-patient batch correction, as recommended). All 49,587 cells received a functional-state label (none below the minimum confidence threshold of 0.2); projection confidence (`functional.cluster.conf`) median 0.814.
- **Cross-validation and identity correction**: crossing `functional.cluster` (ProjecTILs) against `annotation_final` (CellTypist) revealed that 5,293 cells labeled "Regulatory T cells" by CellTypist were projected predominantly to exhausted CD8 states (CD8.TEX 76.6%, CD8.TPEX 11.3%) by ProjecTILs. Raw marker expression on these specific cells (CD8A mean 2.09, 97.1% positive; CD8B mean 1.05, 86.7% positive; vs. CD4 mean 0.097 and FOXP3 mean 0.035, both near-absent despite high nominal detection rates attributable to ambient RNA/dropout noise) confirmed these are genuine CD8 T cells mislabeled by CellTypist, not true regulatory T cells. These cells were reintegrated into the CD8 lineage for downstream analysis. A minor, expected overlap between CD8.TEMRA and NK cells (CD16+/-) was also observed (76.1% of 159 NK cells), consistent with known shared NK-receptor expression in terminally differentiated CD8 T cells, and left uncorrected given the small cell count involved.
- **CD8 composition by pathological response**: within the revised CD8 lineage (42,349 cells, including the reclassified regulatory-T-cell-labeled cells), the proportion of each ProjecTILs functional state was compared between MPR and non-MPR tumors using a global chi-squared test on the full composition table, followed by targeted `prop.test` per individual state; the same statistical approach used for the Bloc5 TME composition analysis. Color coding reused the official ProjecTILs atlas palette for visual and terminological consistency with the reference atlas. p-value display in all barplots (Bloc5 and Bloc6) was standardized using `scales::pvalue()` rather than a manually defined formatting rule.

**LUSC**

- **CD8 purification strategy**: identical method to LUAD (scGate `TME_HiRes`, Harmony reduction supplied to avoid per-level PCA recomputation), applied to the full LUSC TME object (478,223 cells). Runtime ~3h04, longer than the proportional scaling from LUAD would suggest (LUSC is ~1.78x larger but took ~2.4x longer), not indicative of a methodological issue.
- **scGate_multi distribution**: Bcell 87,409; CD4T 135,429; CD8T 108,410; Macrophage 12,631; Mast 9,514; Monocyte 21,649; NK 28,643; panDC 11,064; PlasmaCell 11,827; Neutrophils 7,014; Epithelial 9; Multi (ambiguous) 5,615; unassigned 39,009. The 9 epithelial cells (vs. 0 in LUAD) are consistent with a residual contamination artifact already documented in Bloc5 (cluster 19, Harmony integration artifact), not a concern.
- **Cross-check against CellTypist**: 66.4% agreement (74,133 of 111,596 `lineage_CD8` cells confirmed by scGate), lower than LUAD (72.3%), a histology-level difference noted but not yet interpreted.
- **Cycling/IFN-high cell removal (UCell)**: identical signatures and thresholds as LUAD. 12,177 cells (11.2%, lower than LUAD's 14.8%) removed from 108,410 scGate-confirmed CD8T cells; IFN criterion accounted for ~71% of cumulative flags (vs. ~81% in LUAD), with proportionally more cell-cycle-driven removal in LUSC. 96,233 cells retained for ProjecTILs projection.
- **ProjecTILs functional-state projection**: identical parameters to LUAD (`filter.cells = FALSE`, `split.by = "sampleID"`). All 96,233 cells received a functional-state label; projection confidence median 0.858.
- **Cross-validation and identity correction**: 6,164 cells labeled "Regulatory T cells" by CellTypist were projected predominantly to exhausted CD8 states by ProjecTILs (CD8.TEX 75.5%, CD8.TPEX 7.4%). Raw marker expression on these cells (CD8A mean 2.09, 96.5% positive; CD8B mean 1.17, 81.6% positive; vs. CD4 mean 0.096 and FOXP3 mean 0.034, both near-absent) confirmed the same conclusion as LUAD: genuine CD8 T cells mislabeled by CellTypist. These cells were reintegrated into the CD8 lineage. A smaller, expected overlap between CD8.TEMRA and NK cells was also observed (63.1% of 731 NK cells, vs. 76.1% in LUAD), left uncorrected given the small cell count.
- **CD8 composition by pathological response**: within the revised CD8 lineage (74,357 cells), functional-state composition was compared between MPR and non-MPR using the same statistical approach as LUAD (global chi-squared test, targeted `prop.test` per state, official ProjecTILs atlas color palette, standardized `scales::pvalue()` display).

---

## Preliminary Observations

### Bloc5_Annotation_TME

**LUAD**

- CellTypist (`Immune_All_Low`) directly confirmed 12/17 clusters (≥95% single-family concordance); the remaining 5 clusters were resolved at the cell level. 86.5% of all cells received a marker-confirmed identity; 13.0% remained explicitly labeled "Unresolved" rather than assigned an unverified identity, and 0.6% were excluded as non-immune (lung epithelial) contamination — a pattern independently corroborated by the CellTypist authors' own documentation of an equivalent artifact.
- TME composition differed significantly between MPR and non-MPR tumors (chi-squared test, χ² = 15,402, df = 17, p < 2.2 × 10⁻¹⁶).
- Memory B cells were markedly enriched in MPR tumors (41.2% vs 18.0%), consistent with the presence of a distinct germinal-center B cell population.
- Macrophages were nearly three-fold more abundant in non-MPR tumors (8.1% vs 3.0%; 13.5% vs 6.2% when classical monocytes are included), the largest relative effect size among tested categories.
- Tem/Trm cytotoxic T cells showed a statistically significant but modest difference (20.9% vs 18.8%); across all marker-confirmed CD8 states combined, non-MPR tumors showed a slightly higher overall CD8 proportion (21.9% vs 19.9%).
- 57,782 cells were confirmed as CD8 lineage across all states; 19,238 cells were confirmed as Macrophage lineage (32,450 including classical monocytes, relevant to the TRM/MoMac distinction to be addressed in Bloc7).
- The biological pattern (B cell enrichment in MPR; Macrophage/TAM and cytotoxic CD8 enrichment in non-MPR) is consistent with prior pipeline iterations of this dataset, supporting robustness of the signal independent of the specific annotation method used.

**LUSC**

- CellTypist directly confirmed 13/21 clusters (≥95% single-family concordance); the remaining 8 were resolved at the cell level. 0 cells were left unannotated.
- TME composition differed significantly between MPR and non-MPR tumors (chi-squared test, χ² = 11,861, df = 22, p < 2.2 × 10⁻¹⁶).
- Memory B cells were modestly enriched in MPR tumors (19.1% vs 17.6%), a much smaller effect than the corresponding B-cell enrichment observed in LUAD.
- Tem/Trm cytotoxic T cells showed no detectable difference between MPR and non-MPR (21.3% vs 21.3%, p = 0.92); across all marker-confirmed CD8 states combined, a statistically significant but minimal difference was observed (23.7% vs 23.1%).
- Macrophage/TAM lineage (dominated by Alveolar macrophages) was more abundant in MPR tumors (4.0% vs 3.3% for confirmed Macrophage/TAM alone; 11.6% vs 9.6% including classical monocytes) — the **opposite direction** from LUAD, where Macrophage/TAM was markedly higher in non-MPR tumors.
- 111,596 cells were confirmed as CD8 lineage across all states; 16,920 cells were confirmed as Macrophage/TAM lineage (49,518 including classical monocytes).
- One resolved cluster (391 cells) originated entirely from a single patient and is documented as a Harmony integration artifact rather than a distinct biological subpopulation (see Methods).

**Cross-histology comparison**: the direction of the TAM association with response reverses between LUAD (higher in non-MPR) and LUSC (higher in MPR), while the CD8 association is significant but modest in LUAD and essentially absent in LUSC. This pattern is consistent with the hypothesis that CD8 composition is more conserved across histologies than TAM composition, which appears comparatively more variable, a distinction addressed further in Bloc7 (see Discussion).

### Bloc6_CD8

**LUAD**

- Independent scGate-based CD8 isolation (58,203 cells) showed 72.3% agreement with the CellTypist-derived `lineage_CD8` population (57,782 cells; 41,752 cells confirmed by both methods). This level of concordance is considered acceptable at this stage given the methodological differences between the two approaches, but is not yet interpreted biologically; ProjecTILs projection (CD8-specific reference atlas) is expected to help resolve which discordant cells represent genuine CD8 states versus annotation noise from either method.
- UCell-based removal of cycling and IFN-high cells (official atlas construction thresholds) removed 8,616 cells (14.8%) from the scGate-confirmed CD8T population, leaving 49,587 cells. The IFN criterion accounted for the large majority of flagged cells (~81% of cumulative flags) relative to the cell-cycle criteria — plausibly reflecting an active immune response under neoadjuvant anti-PD-1 treatment, though this is not yet interpreted biologically at this stage.
- ProjecTILs-based functional states were dominated by CD8.EM (32.2%), CD8.TEX (30.9%), and CD8.CM (26.5%), with CD8.MAIT, CD8.TPEX, CD8.TEMRA, and CD8.NaiveLike each below 3.5%. A cross-validation step identified and corrected a CellTypist labeling error affecting 5,293 cells (10.7% of the projected CD8 population) originally annotated as regulatory T cells but confirmed as genuine CD8 T cells (predominantly exhausted/TEX) by both independent projection and raw marker expression (CD8A/CD8B high, CD4/FOXP3 near-absent). This correction is particularly relevant given the central role of the TEX/TPEX axis in this study's neoadjuvant anti-PD-1 context.
- CD8 functional-state composition differed significantly between MPR and non-MPR tumors (chi-squared test, χ² = 259.87, df = 6, p < 0.001). The two largest and most significant effects were CD8.CM, enriched in MPR (6.8 percentage points higher), and CD8.TEX, enriched in non-MPR (6.3 percentage points higher) — consistent with a shift toward a less-differentiated, less-exhausted CD8 phenotype in tumors achieving a major pathological response. CD8.TEMRA, CD8.EM, and CD8.MAIT also reached statistical significance but with substantially smaller effect sizes (under 2 percentage points each), a pattern likely reflecting the statistical power of the large cell numbers involved (~42,000 cells) rather than a comparably strong biological effect; this distinction between statistical and practical significance should be addressed explicitly in the Discussion. CD8.NaiveLike and CD8.TPEX showed no detectable difference between conditions.

**LUSC**

- Independent scGate-based CD8 isolation (108,410 cells) showed 66.4% agreement with the CellTypist-derived `lineage_CD8` population (111,596 cells; 74,133 cells confirmed by both methods) — lower than LUAD (72.3%), consistent with a higher rate of CellTypist "Unresolved" cells observed in LUSC at the cross-validation step (15.2% vs. 9.2% in LUAD).
- UCell-based removal of cycling and IFN-high cells removed 12,177 cells (11.2%, lower than LUAD's 14.8%) from the scGate-confirmed CD8T population, leaving 96,233 cells. The IFN criterion still dominated (~71% of cumulative flags) but proportionally less than in LUAD (~81%), with proportionally more cell-cycle-driven removal in LUSC.
- ProjecTILs-based functional states were dominated by CD8.CM (33.8%), CD8.EM (29.2%), and CD8.TEX (23.4%), a baseline profile shifted toward less-exhausted, more terminally differentiated states compared to LUAD (CD8.EM 32.2%, CD8.TEX 30.9%, CD8.CM 26.5%). A cross-validation step identified and corrected the same CellTypist labeling error as in LUAD, affecting 6,164 cells (8.3% of the projected CD8 population) originally annotated as regulatory T cells but confirmed as genuine CD8 T cells (predominantly exhausted/TEX) by both independent projection and raw marker expression, closely mirroring the LUAD result.
- CD8 functional-state composition differed significantly between MPR and non-MPR tumors (chi-squared test, χ² = 1684.9, df = 6, p < 0.001). The two largest and most significant effects were CD8.TEX, enriched in non-MPR (12.7 percentage points higher), and CD8.CM, enriched in MPR (6.4 percentage points higher), consistent with a shift toward a less-differentiated, less-exhausted CD8 phenotype in tumors achieving a major pathological response. CD8.TEMRA and CD8.EM also reached statistical significance with moderate effect sizes (2 to 3 percentage points). CD8.MAIT, CD8.TPEX, and CD8.NaiveLike reached statistical significance despite substantially smaller effect sizes (under 1.5 percentage points each), likely reflecting the statistical power afforded by the large number of cells analyzed (approximately 74,000) rather than a comparably strong biological effect.
  
**Cross-histology comparison**: comparing CD8 composition changes between MPR and non-MPR within each histology, CD8.CM and CD8.TEX showed the same direction of effect in both LUAD and LUSC (CD8.CM enriched in MPR, CD8.TEX enriched in non-MPR), with the CD8.TEX effect nearly twice as large in LUSC (+12.7 percentage points) as in LUAD (+6.3 percentage points). This is a notably more consistent pattern than the TAM composition result (Bloc5), which reversed direction entirely between the two histologies. CCD8.TPEX itself did not follow this cross-histology consistency: it showed no significant difference in LUAD, but was significantly enriched in non-MPR tumors in LUSC (+0.5 percentage points), the opposite direction from the non-significant trend observed in LUAD. CD8.EM and CD8.MAIT also showed opposite directions of effect between histologies, with small magnitudes in both.

**Going forward**: despite CD8.CM showing an equally strong and consistent compositional shift in both histologies, subsequent functional analyses (UCell scoring, scRepertoire clonality, CollecTRI regulon activity) will focus specifically on the CD8.TEX/CD8.TPEX axis, based on its established mechanistic role in anti-PD-1 response (literature justification to be detailed in the preprint Discussion) rather than on the strength of the compositional signal alone.

---

## Discussion

*To be written once results are available. Condensed prose, distinct from block-level preliminary observations.*
