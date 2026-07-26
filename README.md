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
| Bloc6_CD8 | ProjecTILs CD8 subtyping, UCell signature scoring, CollecTRI regulon activity | To do |
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

- **CD8 purification strategy**: rather than subsetting on the CellTypist-derived `lineage_CD8` flag, CD8 T cells were re-isolated independently on the full TME object of each histology using scGate (`TME_HiRes` model, human), following the same isolation-then-projection strategy used in the official construction of the human CD8 TIL reference atlas (Andreatta et al., 2021, *Nat Commun*; atlas construction pipeline: scGate → UCell cycling/ISG+ removal → STACAS). This avoids circularity with the project's own CellTypist-based annotation and keeps downstream steps methodologically consistent with the reference space used for later ProjecTILs projection.
- **Performance note**: scGate recalculates a dimensionality reduction at each hierarchical gating level by default; supplying the pre-computed Harmony embedding (`reduction = "harmony"`) instead avoids this recomputation and was necessary for tractable runtime on datasets of this size.
- **Cross-check against CellTypist**: `scGate_multi == "CD8T"` (not `lineage_CD8`, and not an intersection/union of the two) was retained as the base population carried forward to the next step, to preserve consistency with the reference atlas construction protocol; `lineage_CD8` is retained as an annotation/cross-validation column only. Concordance between the two methods was assessed but not used to filter cells at this stage; biological interpretation of any discordance is deferred until after ProjecTILs projection, which uses a CD8-specific reference atlas.
- **Cycling/IFN-high cell removal (UCell)**: applied to the scGate-confirmed CD8T population, using the identical signatures and thresholds from the official CD8 human TIL atlas construction pipeline (`cellCycle.G1S`, `cellCycle.G2M` from SignatuR, Tirosh et al., 2016, *Science*; a manually-defined 4-gene IFN signature, *ISG15, IFI6, IFI44L, MX1*, reported by the atlas authors to perform better with UCell than the full SignatuR IFN signature). Thresholds: `cellCycle.G1S < 0.1`, `cellCycle.G2M < 0.1`, `IFN < 0.25`.
- **ProjecTILs functional-state projection**: the purified, non-cycling CD8 population was projected onto the official human CD8 TIL reference atlas (Andreatta et al., 2021, *Nat Commun*; atlas construction: `carmonalab/CD8_human_TIL_atlas_construction`), using `Run.ProjecTILs()` with `filter.cells = FALSE` (query already pre-purified, per official guidance for pre-filtered data) and `split.by = "sampleID"` (per-patient batch correction, as recommended).
- **Cross-validation and identity correction**: `functional.cluster` (ProjecTILs) was crossed against `annotation_final` (CellTypist) to assess concordance between the two independent methods. Where a discrepancy affected a substantial number of cells, identity was arbitrated using raw marker expression (mean log-normalized expression and detection rate) for lineage-defining genes, rather than relying on either method's label alone; detection rate alone was treated with caution given ambient RNA/dropout noise in sparse single-cell data, and mean expression amplitude was used as the primary criterion. Corrections were applied to a CD8-lineage flag used for downstream Bloc6/Bloc7 analyses only, without retroactively altering the Bloc5 `annotation_final` labels or TME composition statistics already computed from them.
- **CD8 composition by pathological response**: within the revised CD8 lineage, per-patient proportions of each ProjecTILs functional state were compared between MPR and non-MPR tumors using the Wilcoxon rank-sum test, with Benjamini-Hochberg correction applied across the seven states. Patients with fewer than 50 CD8 cells were excluded from this analysis, since rare states (such as CD8.TPEX, ~1–3% of CD8) have a high probability of a false zero from under-sampling alone below this threshold. This patient-level approach follows Hu et al. (2023, *Genome Medicine*), a comparable neoadjuvant anti-PD-1 NSCLC scRNA-seq study with a similarly small and unbalanced MPR/non-MPR cohort, and supersedes an earlier cell-level chi-squared/`prop.test` approach (see note below). Color coding reused the official ProjecTILs atlas palette; p-value display was standardized using `scales::pvalue()`.
- **Note on statistical unit (pseudoreplication)**: cell-level tests (chi-squared, `prop.test`) treat each cell as an independent observation despite within-patient correlation, which can inflate statistical power. This was identified while evaluating `propeller` (Phipson et al., 2022, *Bioinformatics*) on the Bloc5 TME composition data, where a cell-level-significant result did not survive patient-level testing. All composition and functional-scoring comparisons in this project were subsequently re-run at the patient level (see above); the Bloc5 TME composition analysis itself was not retested, as it is secondary to the project's central focus on CD8.TEX biology and TAM complexity (Bloc7). A formal interaction test (histology × response) was considered but not pursued, given the absence of precedent in comparable published studies, including the original dataset publication (Liu et al., 2025, *Cell*); LUAD/LUSC comparisons are presented descriptively.
- **Functional scoring within CD8.TEX (UCell)**: to quantify the intensity of specific functional programs beyond the discrete ProjecTILs states, three short published gene signatures (cytotoxicity: *GZMB, PRF1, FASLG*; exhaustion: *HAVCR2, CTLA4, PDCD1, TIGIT, LAG3*; stemness: *TCF7, SELL, LEF1, IL7R*; Carmona et al., 2020, *Oncoimmunology*, retrieved via `SignatuR$Hs$Programs$Tcell.*`) were scored with UCell within the CD8.TEX population, restricted to the revised CD8 lineage. These signatures were chosen over an alternative source (`signaturesHumanTILs`, ProjecTILs atlas-derived state signatures also available from the Carmona lab) because the latter largely reproduce the discrete states already assigned by ProjecTILs itself, whereas the SignatuR programs provide a continuous, independent readout of functional intensity within those states. Per-patient median scores were compared between MPR and non-MPR using the Wilcoxon rank-sum test with Benjamini-Hochberg correction. Patients contributing fewer than 10 CD8.TEX cells were excluded, following the default filtering used in muscat (Crowell et al., 2020), since a median computed over very few cells is unreliable; this retained 47 LUAD patients (9 MPR, 38 non-MPR) and 93 LUSC patients (29 MPR, 64 non-MPR). Note that this threshold differs from the 50-cell minimum used for the composition analysis because the two apply to different denominators: composition requires enough total CD8 cells to compute a proportion, whereas functional scoring requires enough CD8.TEX cells specifically to compute a reliable median within that state.
- **Exclusion of CD8.TPEX from functional analyses**: applying the same 10-cell minimum, only 4 MPR patients would have qualified for CD8.TPEX in LUAD, compared with 13 in LUSC. Analysing this state in LUSC alone would have introduced an asymmetric analytical framework between histologies, precluding any meaningful cross-histology comparison, so CD8.TPEX was excluded from both — neither tested nor reported descriptively, to avoid combining tested and untested results within a single argument. This focus on the terminally exhausted state is also consistent with its reported clinical relevance: baseline CD8.TEX proportion has been shown to be higher in non-responders than responders prior to immunotherapy and to independently predict checkpoint inhibitor resistance (Fang et al., 2026).
- **Transcription factor activity within CD8.TEX (CollecTRI)**: to complement gene-signature scoring with a network-based measure of regulatory activity, transcription factor (TF) activity was inferred using `decoupleR`'s univariate linear model method (`run_ulm`) with the CollecTRI network (Müller-Dott et al., 2023, *Nucleic Acids Research*; ~1,183 human TFs, downloaded directly from static.omnipathdb.org after `decoupleR::get_collectri()` failed due to an unresolved package issue, reformatted to match the expected input, with directionally-ambiguous interactions, 0.4% of the network, excluded). Rather than testing each TF's cell-level activity across patients (which yielded only 5 significant TFs out of 755 tested in an initial analysis, given limited power at n=9-10 MPR patients), a pseudobulk approach was used: raw counts were aggregated per patient within CD8.TEX, a differential expression analysis was run with DESeq2 (MPR vs non-MPR), and `run_ulm` was applied once to the resulting per-gene differential statistic. This approach, recommended in decoupleR's own documentation for this type of comparison, substantially increased sensitivity (67 of 637 TFs significant in LUAD, 10.5%, versus 14 of 649 in LUSC, 2.2%) without introducing repeated per-TF testing across patients. No TFs were pre-selected; all those passing the network's minimum target-gene coverage (`minsize = 5`) were tested.
- **CD8.TPEX and TF activity**: consistent with its exclusion from signature scoring, CD8.TPEX was not analysed with CollecTRI.

---

## Preliminary Observations

### Bloc5_Annotation_TME

**LUAD**

- CellTypist (`Immune_All_Low`) directly confirmed 12/17 clusters (≥95% single-family concordance); the remaining 5 clusters were resolved at the cell level. 86.5% of all cells received a marker-confirmed identity; 13.0% remained explicitly labeled "Unresolved" rather than assigned an unverified identity, and 0.6% were excluded as non-immune (lung epithelial) contamination — a pattern independently corroborated by the CellTypist authors' own documentation of an equivalent artifact.
- 57,782 cells were confirmed as CD8 lineage across all states; 19,238 cells were confirmed as Macrophage lineage (32,450 including classical monocytes, relevant to the TRM/MoMac distinction to be addressed in Bloc7).

**LUSC**

- CellTypist directly confirmed 13/21 clusters (≥95% single-family concordance); the remaining 8 were resolved at the cell level. 0 cells were left unannotated.
- 111,596 cells were confirmed as CD8 lineage across all states; 16,920 cells were confirmed as Macrophage/TAM lineage (49,518 including classical monocytes).
- One resolved cluster (391 cells) originated entirely from a single patient and is documented as a Harmony integration artifact rather than a distinct biological subpopulation (see Methods).

**Note**: formal MPR vs non-MPR TME composition comparisons (cell-level chi-squared) were computed but are not reported here, as patient-level re-analysis (propeller) did not confirm robust significance at the available sample size (10 MPR vs 42 non-MPR in LUAD), and this axis is secondary to the project's central focus on CD8 TEX/TPEX biology and TAM complexity (Bloc7).

### Bloc6_CD8

**CD8.TPEX rarity (both histologies)**, CD8.TPEX did not receive a formal patient-level significance test, given insufficient cell counts per patient in both histologies (unrelated to pathological response; zero-cell rates were similar between MPR and non-MPR). CD8.TEX, in contrast, supported a full patient-level test.

| | LUAD | LUSC |
|---|---|---|
| CD8.TEX — % patients with ≥1 cell | 98.1% | 95.2% |
| CD8.TEX — of those, % with <10 cells | 7.8% | 7.0% |
| CD8.TPEX — % patients with ≥1 cell | 88.5% | 82.9% |
| CD8.TPEX — of those, % with <10 cells | 43.5% | 60.9% |

No published study comparing CD8.TPEX abundance between LUAD and LUSC was identified.

**LUAD**

- Independent scGate-based CD8 isolation (58,203 cells) showed 72.3% agreement with the CellTypist-derived `lineage_CD8` population (57,782 cells; 41,752 cells confirmed by both methods). This level of concordance is considered acceptable at this stage given the methodological differences between the two approaches, but is not yet interpreted biologically; ProjecTILs projection (CD8-specific reference atlas) is expected to help resolve which discordant cells represent genuine CD8 states versus annotation noise from either method.
- UCell-based removal of cycling and IFN-high cells (official atlas construction thresholds) removed 8,616 cells (14.8%) from the scGate-confirmed CD8T population, leaving 49,587 cells. The IFN criterion accounted for the large majority of flagged cells (~81% of cumulative flags) relative to the cell-cycle criteria — plausibly reflecting an active immune response under neoadjuvant anti-PD-1 treatment, though this is not yet interpreted biologically at this stage.
- ProjecTILs-based functional states were dominated by CD8.EM (32.2%), CD8.TEX (30.9%), and CD8.CM (26.5%), with CD8.MAIT, CD8.TPEX, CD8.TEMRA, and CD8.NaiveLike each below 3.5%. A cross-validation step identified and corrected a CellTypist labeling error affecting 5,293 cells (10.7% of the projected CD8 population) originally annotated as regulatory T cells but confirmed as genuine CD8 T cells (predominantly exhausted/TEX) by both independent projection and raw marker expression (CD8A/CD8B high, CD4/FOXP3 near-absent). This correction is particularly relevant given the central role of the TEX/TPEX axis in this study's neoadjuvant anti-PD-1 context.
- At the patient level (Wilcoxon rank-sum test, BH-corrected across the seven states, patients with fewer than 50 CD8 cells excluded, 49/52 patients retained: 10 MPR, 39 non-MPR), no functional state reached statistical significance for the MPR versus non-MPR comparison (all BH-adjusted p ≥ 0.91). This includes CD8.TEX, which had appeared significant in an earlier cell-level analysis; the patient-level result is considered more reliable given the risk of pseudoreplication in cell-level testing (see Methodology), and the discrepancy likely reflects limited statistical power at this sample size (only 10 MPR patients) rather than the cell-level result being biologically meaningful.
- Within CD8.TEX, UCell-based scoring of cytotoxicity, exhaustion, and stemness gene programs (Carmona et al., 2020) at the patient level revealed that CD8.TEX cells in non-MPR tumors were both more cytotoxic (median 0.40 versus 0.29; BH-adjusted p = 0.004) and more exhausted (median 0.32 versus 0.18; p = 0.001) than those in MPR tumors, a pattern not mutually exclusive given that exhausted CD8 T cells are known to retain partial cytotoxic potential. The stemness score was near zero in both groups and did not differ significantly (p = 0.665), likely reflecting a known technical limitation of short gene signatures with UCell rather than a biological effect. Notably, this functional-scoring result contrasts with the corresponding composition analysis, where no CD8 functional state reached significance in LUAD at the patient level — suggesting pathological response may be more closely associated with the functional intensity of CD8.TEX cells than with their proportion within the CD8 compartment.
- Transcription factor activity inference (CollecTRI, pseudobulk + DESeq2) identified 67 of 637 tested transcription factors with significantly different activity between MPR and non-MPR tumors (10.5 percent). The strongest and most significant associations were with NF-κB (RELA, NFKB1), STAT family members (STAT1, STAT5B, STAT3), AP-1 components (JUN, FOS), HSF1, and EGR1, all more active in non-MPR tumors; only four transcription factors (ZBTB4, ZBTB14, YBX1, ZBED1) showed significantly higher activity in MPR tumors. This asymmetry mirrors the composition result, where CD8.TEX proportion did not differ significantly between conditions, suggesting pathological response is associated with the transcriptional program of TEX cells rather than their relative abundance.

**LUSC**

- Independent scGate-based CD8 isolation (108,410 cells) showed 66.4% agreement with the CellTypist-derived `lineage_CD8` population (111,596 cells; 74,133 cells confirmed by both methods), lower than LUAD (72.3%), consistent with a higher rate of CellTypist "Unresolved" cells observed in LUSC at the cross-validation step (15.2% vs. 9.2% in LUAD).
- UCell-based removal of cycling and IFN-high cells removed 12,177 cells (11.2%, lower than LUAD's 14.8%) from the scGate-confirmed CD8T population, leaving 96,233 cells. The IFN criterion still dominated (~71% of cumulative flags) but proportionally less than in LUAD (~81%), with proportionally more cell-cycle-driven removal in LUSC.
- ProjecTILs-based functional states were dominated by CD8.CM (33.8%), CD8.EM (29.2%), and CD8.TEX (23.4%), a baseline profile shifted toward less-exhausted, more terminally differentiated states compared to LUAD (CD8.EM 32.2%, CD8.TEX 30.9%, CD8.CM 26.5%). A cross-validation step identified and corrected the same CellTypist labeling error as in LUAD, affecting 6,164 cells (8.3% of the projected CD8 population) originally annotated as regulatory T cells but confirmed as genuine CD8 T cells (predominantly exhausted/TEX) by both independent projection and raw marker expression, closely mirroring the LUAD result.
- At the patient level (same approach, 101/105 patients retained: 34 MPR, 67 non-MPR), CD8.TEX proportion was significantly higher in non-MPR tumors (median 25.8% versus 11.2%; BH-adjusted p = 0.0135). CD8.CM and CD8.MAIT proportions were also significantly higher in MPR tumors (BH-adjusted p = 0.034 for both). All other states, including CD8.TPEX, were not significant after correction.
- Within CD8.TEX, UCell-based scoring at the patient level showed the same pattern as LUAD: significantly higher cytotoxicity (median 0.40 versus 0.31; BH-adjusted p < 0.001) and exhaustion (median 0.31 versus 0.22; p = 0.001) in non-MPR tumors, with stemness scores remaining near zero and not differing significantly (p = 0.175). This patient-level result was more consistent between histologies than the corresponding composition analysis, which had reached significance in LUSC but not in LUAD.
- Transcription factor activity inference identified only 14 of 649 tested transcription factors as significant (2.2 percent), despite a larger patient cohort than LUAD (93 versus 47 patients) — a point to revisit in the Discussion once compared against the literature, rather than assumed to reflect a lack of statistical power. The dominant pattern involved E2F family members and TFDP1 (E2F4, E2F2, TFDP1, E2F1, associated with cell-cycle regulation), together with IRF3 and TBP. STAT1 and EGR1 were the only two transcription factors reaching significance in both histologies; SNAI2 was the only transcription factor significantly more active in MPR tumors.
  
**Summary**: across all three analyses of the CD8 lineage (composition, functional scoring, and transcription factor activity), CD8.TEX consistently showed a significant association with pathological response in both histologies for functional scoring and transcription factor activity, while the compositional association was significant in LUSC only. CD8.TPEX was excluded from functional analyses in both histologies due to insufficient cell numbers per patient (see Methods).


---

## Discussion

*To be written once results are available. Condensed prose, distinct from block-level preliminary observations.*
