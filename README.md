# LUAD vs LUSC: CD8 and TAM Biology in Neoadjuvant Anti-PD-1 Treated NSCLC

Preprint: [DOI to be added once online]

## Research Question

How does CD8 and TAM biology (transcriptional activity, subtype heterogeneity, cell-cell communication) differ between lung adenocarcinoma (LUAD) and lung squamous cell carcinoma (LUSC) in patients treated with neoadjuvant anti-PD-1 therapy?

## Dataset

- **Source**: GSE243013 (Liu et al., 2025, *Cell*), single-cell RNA/TCR-seq of surgical tumor samples after neoadjuvant chemo-immunotherapy in NSCLC.

- **Cohort**: the GEO metadata includes 243 sequenced patients. Nine received chemotherapy alone and are excluded at Script 00 (Bloc0), since this study concerns anti-PD-1 response specifically. This leaves 234 patients (61 LUAD, 173 LUSC), verified by unique `sampleID` and cross-checked against the `pathological_response` distribution (MPR n=42, non-MPR n=108, pCR n=83, unknown n=1).

- **Analysed patients**: pCR patients and the single patient with an unclassified response are excluded in the last script of Bloc1, so that all downstream analyses compare MPR with non-MPR. The analysed cohort comprises 150 patients (50 LUAD: 10 MPR, 40 non-MPR; 100 LUSC: 32 MPR, 68 non-MPR). The full patient flow is given in the table of Bloc1_QC.

- **Metadata file**: `GSE243013_NSCLC_immune_scRNA_metadata.csv.gz`, cell-level resolution, key columns include `sampleID`, `cancer_type`, `pathological_response`, `major_cell_type`, `sub_cell_type`.

- **EGFR status note:** Not available in public metadata as a dedicated field. No EGFR TKI recorded in the `targeted_therapy` field for either arm. Patients with EGFR mutation or ALK rearrangement were reported as excluded by design in the source publication (Liu et al., 2025), since they would be eligible for targeted therapies instead, though this could not be verified from the metadata available here.

---

## Summary of findings

Results by block are given further down; this section summarizes what the analyses converged on.

**Abundance of the major populations largely does not separate response groups, but the state of CD8 cells does**. At the patient level, the proportions of macrophages and of cytotoxic Tem/Trm cells (the display-level proxy for CD8 cells) do not differ between responders and non-responders in either histology. The only difference in the composition of the microenvironment is a higher proportion of regulatory T cells in non-responders in squamous carcinoma (median 18.0 versus 9.81 percent; adjusted P = 2.6e-4). Within the exhausted CD8 compartment, CD8.TEX cells of non-responders score higher for both cytotoxicity and exhaustion in both histologies (adjusted P from 0.001 to 0.008), and are more abundant in non-responders in squamous carcinoma only (median 22.6 versus 9.0 percent; adjusted P = 0.005). These are transcriptional scores and do not measure cytotoxic capacity. Clonal expansion and diversity of the TCR repertoire within CD8.TEX did not differ between response groups in either histology.

**Regulatory activity does separate them in most populations, and almost entirely on one side.** In the myeloid compartment, the factors reaching significance were more active in non-responders in both populations and both histologies: in LUAD, 62 of 64 in Macrophage_pure and 33 of 34 in Monocyte_pure; in LUSC, 7 of 8 and 8 of 11. In exhausted CD8 cells, 7 factors differed in LUAD, all more active in non-responders, and none reached significance in LUSC.

**Known axes described in the literature are found in the output of both histologies, mostly on the non-responder side.** Of nine axes described in the literature for CD8 and myeloid cells, all nine appear in the prioritisation output for non-responders in LUAD and eight in LUSC. In responders, only CXCL16–CXCR6 appears in LUAD, and CXCL16–CXCR6 and SPP1–CD44 in LUSC. Presence in the output indicates neither a stronger interaction nor a significant difference.

**Transcription factor inference in exhausted CD8 cells depends on how well the network documents the differentially expressed genes.** In LUAD, 65 genes differed between response groups, of which 24 (37%) have no documented regulator in the network used (CollecTRI), and 7 factors reached significance. In LUSC, 258 genes differed, of which 146 (57%) have no documented regulator, and no factor reached significance. The absence of significant factors in LUSC likely reflects, at least in part, this limited network coverage rather than an absence of regulatory differences between response groups.

**What these results cannot establish.** Sampling is confined to a single post-treatment timepoint, so a microenvironment that contributed to treatment failure cannot be distinguished from one that persisted alongside residual disease. Interactions are inferred from co-expression without spatial information. Patients achieving a complete pathological response were not included in the comparison.

---

## Repository Structure

```
TAM_CD8_LUAD_LUSC_scRNAseq/
├── Data/                                 # input data (not tracked, download from GEO, see Dataset)
├── Objects/                              # saved Seurat/R objects (not tracked, too large; regenerated by running the scripts in order)
├── References/                           # ProjecTILs human CD8 TIL reference atlas (Andreatta et al., 2021), not tracked
├── scGateDB/                             # scGate models downloaded with get_scGateDB() (not tracked)
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
├── Results/
│   ├── Bloc0_DataPrep/{Figures,Tables}
│   ├── Bloc1_QC/{Figures,Tables}
│   ├── Bloc2_Norm_HVG_PCA/{Figures,Tables}
│   ├── Bloc3_Harmony/{Figures,Tables}
│   ├── Bloc4_Clustering_UMAP/{Figures,Tables}
│   ├── Bloc5_Annotation_TME/{Figures,Tables}
│   ├── Bloc6_CD8/{Figures,Tables}
│   ├── Bloc7_TAM/{Figures,Tables}
│   └── Bloc8_Communication/{Figures,Tables}
├── Exploratory_analyses.md               # analyses not included in the preprint
├── renv.lock                             # R package versions (renv)
├── requirements.txt                      # Python packages (CellTypist)
├── LICENSE                               # MIT
└── README.md
```

**Naming convention**: within each `Scripts/BlocX_*/` folder, scripts are numbered starting at `01`, following the execution order *within that block* (not a global pipeline count). Most scripts carry a histology suffix, `LUAD` or `LUSC`, and are run separately on each histology. Example: `01_LUAD_QC.R` and `01_LUSC_QC.R` in `Bloc1_QC`. Scripts that process both histologies in a single run carry no suffix.

Results (figures and tables) mirror this structure exactly: same block folder, same script number, same histology suffix, so that a given script and its output can always be matched by name alone.

All paths are relative, managed via the `here` package (`here()` / `here("Data", ...)`), for portability across Windows, WSL, and Docker without modification.

---

## Pipeline

| Block | Content | Status |
| --- | --- | --- |
| Bloc0_DataPrep | Matrix loading (BPCells), LUAD/LUSC separation, patient-level clinical metadata | Done |
| Bloc1_QC | Quality control plots and filtering (per-arm thresholds) | Done |
| Bloc2_Norm_HVG_PCA | Normalization, highly variable genes, PCA | Done |
| Bloc3_Harmony | Batch integration | Done |
| Bloc4_Clustering_UMAP | Clustering, UMAP visualization | Done |
| Bloc5_Annotation_TME | TME annotation with CellTypist | Done |
| Bloc6_CD8 | ProjecTILs CD8 subtyping, UCell signature scoring, CollecTRI regulon activity, TCR repertoire | Done |
| Bloc7_TAM | TAM purification (scGate), composition of Macrophage_pure and Monocyte_pure, CollecTRI regulon activity; exploratory scoring analyses are described in `Exploratory_analyses.md` | Done |
| Bloc8_Communication | MultiNicheNet differential cell-cell communication | Done |

Each block's methodology, chosen parameters and rationale are documented below.

## Reproducibility

- **R environment**: package versions are recorded with `renv` (`renv.lock`). Restore them with `renv::restore()`.
- **Python environment**: the CellTypist annotation step of Bloc5 uses Python. Package versions are listed in `requirements.txt` (for example `celltypist==1.7.1`, `scanpy==1.12.2`). Install them with `pip install -r requirements.txt`.
- **Data**: raw data are not included in the repository because of their size. They are obtained from GEO (GSE243013, see Dataset above); Script 00 of Bloc0 builds the reference ID lists used by all downstream scripts.
- **Intermediate objects**: saved objects (`Objects/`) are not versioned either. They are regenerated by running the scripts in order, block by block and by script number (see Repository Structure).
- **Memory**: the analyses were run with 32 GB of RAM. The scGate step on the LUSC object exceeded this and was therefore run in batches of patients (see Bloc6).
- **Randomness**: seeds are fixed in the scripts wherever a random step is involved (for example the assignment of patients to scGate batches in Bloc6).

---

## Tools and Methods

| Tool | Role in the project | Reference | Block |
| --- | --- | --- | --- |
| BPCells | Disk-backed matrix loading, avoids RAM overload | Parks and Greenleaf, 2025 *bioRxiv* | Bloc0 |
| Seurat | Single-cell data handling | Hao et al., 2024, *Nature Biotechnology* | Bloc1-Bloc8 |
| Harmony | Batch integration | Korsunsky et al., 2019, *Nature Methods* | Bloc3 |
| CellTypist | TME annotation (cell-level classification, majority voting) | Domínguez Conde et al., 2022, *Science* | Bloc5 |
| scGate | Marker-based purification of CD8 and macrophage compartments | Andreatta et al., 2022, *Bioinformatics* | Bloc6, Bloc7 |
| ProjecTILs | CD8 subtyping via reference atlas projection | Andreatta et al., 2021, *Nature Communications* | Bloc6 |
| UCell | Gene signature scoring (CD8: cycling and IFN-high cell removal, functional scoring within CD8.TEX) | Andreatta & Carmona, 2021, *Computational and Structural Biotechnology Journal* | Bloc6 |
| decoupleR + CollecTRI | Transcription factor activity inference (CD8 and myeloid) | decoupleR: Badia-i-Mompel et al., 2022, *Bioinformatics Advances*; CollecTRI: Müller-Dott et al., 2023, *Nucleic Acids Research* | Bloc6, Bloc7 |
| edgeR | Gene filtering before differential expression (`filterByExpr`) | Chen et al., 2016, *F1000Research* | Bloc6, Bloc7 |
| DESeq2 | Pseudobulk differential expression | Love et al., 2014, *Genome Biology* | Bloc6, Bloc7 |
| MultiNicheNet | Differential cell-cell communication, multi-sample multi-condition design | Browaeys et al., 2023, *bioRxiv*; builds on NicheNet: Browaeys et al., 2020, *Nature Methods* | Bloc8 |
| TCR repertoire analysis | Clonal expansion and diversity within CD8.TEX; clonotypes from the authors' annotation file, metrics computed explicitly | Expansion threshold: Liu et al., 2025, *Cell*; diversity indices: Mika et al., 2025 | Bloc6 |
| muscat | Pseudobulk framework; source of the default 10-cell minimum per patient and sample | Crowell et al., 2020, *Nature Communications* | Bloc6, Bloc7, Bloc8 |

---

## Methodology

### Bloc0_DataPrep

- **Exclusion of patients without immunotherapy**: inspection of the clinical metadata (`anti-PD1_therapy` column) identified 9 patients who received chemotherapy alone (`anti-PD1_therapy == "No"`: 2 LUAD, 7 LUSC), confirmed by a cross-check in the console. This is consistent with the numbers reported by Liu et al. 2025 (343 patients followed, 243 with scRNA/TCR-seq: 234 anti-PD-1 + chemotherapy, 9 chemotherapy alone). These 9 patients are excluded at Script 00, before the LUAD/LUSC split; the script halts if the number of such patients is not exactly 9, and the excluded patients are logged (`00_excluded_chemo_only_patients.csv`). The corrected cohort comprises 61 LUAD and 173 LUSC patients (234 in total).

- **Separation**: LUAD/LUSC split performed once (Script 00), producing reference ID lists (`sampleID`, `cellID`) reused by all downstream scripts. No re-filtering of `cancer_type` elsewhere in the pipeline.

- **Verification**: patient counts checked against previously validated values (61 LUAD, 173 LUSC) at each step; script halts on mismatch.

- **Matrix loading**: BPCells used for disk-backed loading, avoiding full in-memory load of the raw count matrix.

- **Clinical metadata**: one row per patient, deduplicated by `sampleID`. A known typo (`"unknowm"` instead of `"unknown"`) in `pathological_response` was kept as-is for traceability rather than silently corrected.

- **Response categories**: LUAD has 3 response categories (MPR, non-MPR, pCR). LUSC has 4, including one patient labeled `"unknowm"` (unclassified). At Bloc0, all patients of the corrected cohort are kept in the reference ID lists, including this patient and the pCR patients. They are removed at the response group filtering step (last script of Bloc1, before normalisation), so this patient does not enter any downstream score calculation or MPR vs non-MPR comparison.

- **Total cell count cross-check**: the full metadata table contains 1,254,749 cells, of which 41,612 belong to the 9 excluded chemotherapy-only patients. The remaining 1,213,137 cells match LUAD (328,528) + LUSC (884,609).

### Bloc1_QC

- **Metrics**: nFeature_RNA, nCount_RNA, and %.mt (mitochondrial content, `^MT-` pattern) computed per cell, after minimal Seurat object creation filters (`min.cells=3`, `min.features=200`).

- **Threshold selection**: thresholds based on the 99th %ile of each metric's distribution, computed separately for LUAD and LUSC (never copied across histologies), and uniform across all patients and pathological response groups within each histology. Differentiating thresholds by response group was explicitly avoided, since `pathological_response` is the primary comparison axis (MPR vs non-MPR); group-specific thresholds risk introducing an artificial bias into that comparison.

- **Thresholds applied**:
  - LUAD: nFeature_RNA [200-5000], nCount_RNA [1000-28000], %.mt < 10%
  - LUSC: nFeature_RNA [200-4300], nCount_RNA [1000-26000], %.mt < 10%

- **Prior filtering by the source authors**: the cells of GSE243013 had already been filtered by Liu et al. (2025) before deposit, retaining cells with 1,200 to 40,000 UMI counts, 600 to 6,000 detected genes and less than 15% mitochondrial counts, after removal of doublets with Scrublet. As a consequence, the Seurat object creation filters (`min.cells=3`, `min.features=200`) and the lower thresholds applied here (nFeature_RNA > 200, nCount_RNA > 1000) removed no cell: 328,528 cells in LUAD and 884,609 in LUSC before filtering, identical to the loaded matrices. The cells removed (2.3% in LUAD, 2.8% in LUSC) come from the upper thresholds (99th percentile) and from the stricter mitochondrial threshold (10% against 15%). The lower bound of 1,000 UMI counts was chosen without a published reference.

- **Verification**: cells removed per pathological response group checked post-filtering to confirm no group was disproportionately affected. LUAD: MPR -1.8%, non-MPR -2.7%, pCR -1.4%. LUSC: MPR -3.1%, non-MPR -3.0%, pCR -2.5%, unclassified ("unknowm") -2.5%. Comparable removal rates across groups in both histologies support the validity of uniform thresholds for the downstream MPR vs non-MPR comparison.

- **Object portability**: raw counts converted from BPCells (disk-backed) to standard in-memory sparse matrix (`dgCMatrix`) at the end of QC filtering, ensuring downstream objects remain portable across operating systems/environments, independent of the original BPCells storage path.

- **Response group filtering**: following QC filtering, cells are restricted to `pathological_response %in% c("MPR", "non-MPR")`. pCR patients (and the single `"unknowm"`-labeled LUSC patient) are excluded from this point onward in both histologies, since neither is used in any downstream scoring or comparison; the proportion of pCR patients per histology is reported separately, computed from the clinical metadata rather than from the filtered Seurat object. Cell counts after filtering: LUAD 260,674 cells (MPR 50,914 / non-MPR 209,760, from 320,913 before filtering); LUSC 457,403 cells (MPR 170,668 / non-MPR 286,735, from 859,903 before filtering).

- **Patient numbers from the corrected cohort to the analysed object**: the corrected cohort (234 patients) excludes the 9 patients who received no immunotherapy (2 in LUAD, 7 in LUSC), removed at Bloc0 (Script 00). pCR patients (and the single `"unknowm"`-labeled LUSC patient) are then excluded in the last script of Bloc1, before normalisation, so that everything downstream compares MPR with non-MPR patients only.

| | LUAD | LUSC | Total |
|---|---|---|---|
| Corrected cohort (Bloc0, after removal of the 9 patients without immunotherapy) | 61 | 173 | 234 |
| pCR patients (excluded, Bloc1) | 11 | 72 | 83 |
| `"unknowm"` response (excluded, Bloc1) | 0 | 1 | 1 |
| **Patients analysed (MPR + non-MPR)** | **50** | **100** | **150** |
| of which MPR | 10 | 32 | 42 |
| of which non-MPR | 40 | 68 | 108 |

### Bloc2_Norm_HVG_PCA

- **Normalization**: LogNormalize (scale factor 10,000), standard method.

- **Highly variable genes**: 3,000 HVGs for both LUAD and LUSC (vst method), above the Seurat default of 2,000. Verified against `major_cell_type` metadata (CD45+ sorted dataset, no epithelial compartment), intended to better capture fine-grained heterogeneity within immune subpopulations. LUSC initially used a reduced HVG count (2,000) as a temporary computational constraint on the larger pre-filtering dataset; unified to 3,000 (matching LUAD) once the response-group filtering (see Bloc1_QC) reduced the LUSC cell count enough to remove this constraint.

- **Cell cycle regression**: not applied to either histology. Initially applied to LUAD only, then reassessed: regression reduced but did not fully resolve a mixed-lineage proliferation cluster (a residual mixed cluster still required post-hoc resolution via SingleR - see Bloc5), so it did not provide a complete solution even where it succeeded technically. On LUSC, regression additionally caused repeated memory allocation failures on the larger dataset. Since post-hoc SingleR resolution alone proved sufficient on LUAD, regression was dropped entirely for consistency across both histologies.

- **Scaling**: applied to HVGs only, not the full gene set, to limit memory usage. The dense scaled matrix (`scale.data`) is discarded from the object immediately after PCA, since it is not required by any downstream step (Harmony, clustering, UMAP, `FindAllMarkers`, SingleR all operate on the PCA/Harmony embedding or on the normalized sparse matrix).

- **PCA**: 50 components computed for both arms; number of PCs to retain determined from per-PC variance explained, using the point where successive PCs contribute a near-negligible marginal variance (<0.1%).

- **PC1 diagnostic**: correlation with nCount_RNA checked to rule out a purely technical origin for the dominant component (LUAD: r=-0.50; LUSC: r=-0.49 - both moderate, not indicative of a purely technical artifact); correlation with %.mt negligible in both arms (LUAD r≈-0.02; LUSC r≈0.07).

- **PCs retained**: 20 for LUAD (variance plateau from ~PC12-13, with a minor rebound around PC14-16 before flattening again through PC19-20; cumulative variance 81.9% at PC20) and 20 for LUSC (variance plateau from ~PC11-12, with a minor rebound around PC13-14 before flattening again through PC17-20; cumulative variance 82.4% at PC20), slightly beyond the irregular plateau to reduce the risk of losing signal from rare subpopulations of interest (TAM, CD8 substates)

### Bloc3_Harmony

- **Batch correction**: Harmony (`group.by.vars = "sampleID"`), correcting for inter-patient technical variation, applied on the 20 PCs retained in Bloc2.

- **Verification**: patient-level dispersion (SD of per-patient mean coordinate) compared before/after Harmony as a proxy for integration quality. LUAD: 5.95 to 5.43 (-8.7%). LUSC: 5.60 to 5.18 (-7.5%). Comparable reductions across both histologies support the robustness of the integration step.

### Bloc4_Clustering_UMAP

- **Clustering**: `FindNeighbors`/`FindClusters` on the Harmony embedding (not raw PCA), 20 dimensions, resolution 0.5 as a conservative starting point for global TME clustering (20 clusters in both LUAD and LUSC).

- **UMAP**: computed on the Harmony embedding, same 20 dimensions.

- **Integration quality check**: patient distribution visually inspected on the UMAP; well-mixed colors across all major cell islands in both LUAD and LUSC.

### Bloc5_Annotation_TME

*Script order (Scripts/Bloc5_Annotation_TME/): 01 exploratory manual annotation (FindAllMarkers, canonical markers, sanity-check only, not merged); 02 CellTypist export/import (02_LUAD_SingleR.R deprecated, traceability only); 03 mixed-cluster resolution, SingleR (03_LUAD_MixedClusters_Resolution_SingleR.R deprecated, traceability only); 04 mixed-cluster resolution, CellTypist (per-label marker verification, all 20 clusters both histologies, supersedes 95%-homogeneity threshold); 05 final annotation build (annotation_final, lineage_CD8, lineage_TAM, lineage_TAM_related); 06 TME visualization and patient-level testing (Wilcoxon, Benjamini-Hochberg across every confirmed population); 07 patient-level composition of the TAM and CD8 populations (figure panel).*

**Methodology (common to both histologies):**

- **Marker-based annotation**: `FindAllMarkers` (Wilcoxon test via `presto`), top10 markers per cluster reviewed manually against canonical immune markers. 20 clusters (filtered MPR/non-MPR object, no cell-cycle regression) reviewed as an initial exploratory pass; not merged into the final annotation (see below).

- **Cell-level annotation (CellTypist)**: `Immune_All_Low` model (Domínguez Conde et al., 2022, *Science*; multi-tissue atlas including tumor-infiltrating lymphocytes), run with `majority_voting = True`, producing a `CellTypist_majority` label per cell. Applied to all 20 clusters per histology.

- **Minimum group size for marker verification**: within each cluster, `CellTypist_majority` subgroups below 30 cells were excluded from the `FindAllMarkers` comparison — below this size, a Wilcoxon test (the default test underlying `FindAllMarkers`) lacks sufficient statistical power to yield a reliable result. When a cluster retains only one subgroup above this threshold, no comparison group is available and `FindAllMarkers` cannot be run; marker verification then relies on direct expression fetching (`FetchData`), computing the proportion of cells expressing each official marker directly.

- **Cluster-level marker verification**: to avoid the subjectivity of choosing a grouping level for a percentage-based threshold (e.g. CD4 vs CD8, Treg vs Th, macrophage vs monocyte vs DC), every `CellTypist_majority` subgroup was verified individually, in every cluster of both histologies, regardless of its proportion within the parent cluster. Cells were grouped by their individual `CellTypist_majority` label (labels below a minimum-cell threshold excluded as noise), then a targeted `FindAllMarkers` was run per subgroup. Each subgroup's marker profile was checked against the CellTypist authors' own curated marker sets (Mullan et al., 2023, *Front Immunol* for T-cell subtypes; Domínguez Conde et al., 2022 for NK/myeloid/B) — a subgroup was confirmed only if its official marker(s) were recovered among its own top differentially expressed genes; unconfirmed subgroups were retained as explicitly traceable "Unresolved" labels rather than assigned an unverified identity. Signals dominated by non-immune markers (lung epithelial contamination: *SCGB1A1*, *SFTPB*, *SFTPC*, *SCGB3A1*) were excluded as non-immune contamination — a pattern independently reported by the CellTypist authors as an ambient-RNA/doublet artifact in their own macrophage cluster analysis.

- **Final annotation** (`annotation_final`): built entirely from CellTypist. Each retained label is the model's own terminology, not translated into project-specific nomenclature, to preserve direct traceability to the source classifier. Every cell carries a label. The CD8 lineage flag (`lineage_CD8`) was built from confirmed labels only; the TAM and monocyte flags (`lineage_TAM`, `lineage_TAM_related`) were built from the raw CellTypist label regardless of marker-confirmation status, so that transitional or non-canonical states remain available under these flags. These technical flags are used only to seed the dedicated CD8 (Bloc6) and TAM (Bloc7) extractions, not for the statistical comparison described below.

- **Display simplification**: for figure readability only, categories below 500 cells are grouped into a single "Minor/unresolved categories" bin, and non-immune contamination into a separate bin (kept visible regardless of size); this applies to the summary figure (`annotation_display`) only and does not affect `annotation_final`, which retains every individually confirmed or unconfirmed label.

- **Statistical comparison (patient-level, whole TME)**: per-patient proportions of every displayed category (`annotation_display`, excluding the contamination and minor/unresolved bins) were compared between MPR and non-MPR using the Wilcoxon rank-sum test, with Benjamini-Hochberg correction applied across all categories tested together — not restricted in advance to any subset, to avoid the selection bias of choosing which populations to test after seeing the data. Patients contributing fewer than 50 annotated cells were excluded: a proportion computed over a very small denominator is unstable regardless of the population's abundance, since one cell out of twelve already represents eight %. This supersedes an earlier cell-level approach (global chi-squared plus targeted `prop.test`), which reported significant differences that did not survive patient-level testing — see the note on pseudoreplication below. CD8 and TAM are not singled out at this stage: each is retested in its own dedicated block (Bloc6, Bloc7) on its more rigorously defined population (post-scGate/ProjecTILs), where the corresponding figure panel is produced; a composite figure combining both is assembled once both blocks are complete.

- **Composition figure**: one panel per histology — the cumulative composition barplot, retained as a descriptive overview with no test attached to it directly, since cells from one patient are not independent observations; the accompanying patient-level test (above) is reported as a table, not overlaid on this panel. Cell-level testing would also let patients with disproportionately many cells of a given type dominate the result irrespective of how many patients actually support the trend — an uneven contribution visible directly in this dataset, where cluster size per patient varies by two orders of magnitude for the same cell state. Per-patient boxplots for CD8 T cells and macrophages are produced separately in their dedicated blocks (Bloc6, Bloc7) and assembled into a composite figure at that stage.

**Note:** CellTypist was selected for the Bloc5 TME annotation because independent literature validation was readily available: the CellTypist authors' reference profiles cover most major TME lineages (including NK cells), and Mullan et al. (2023) provided a complementary, T-cell-specific reference. This enabled marker-based confirmation of cluster identity (FindAllMarkers against both references) rather than relying on classifier output alone; clusters lacking marker support in either reference were labelled Unresolved. scGate (used in Bloc6-7 for CD8 and TAM purification) is explicitly designed to operate without reference gene expression profiles (Andreatta et al., 2022) — a different validation logic, better suited to targeted purification of a single population than to a first-pass annotation of a heterogeneous tissue.

### Bloc6_CD8

*Script order (`Scripts/Bloc6_CD8/`): 01 scGate CD8 purification, batched (01a whole-object reconstruction, LUAD only; 01b batch verification, full cross-check LUAD only, structural checks both histologies); 02 UCell cycling/IFN-high removal; 03 ProjecTILs state projection (03b per-patient state testability); 04 deprecated (cross-validation vs CellTypist, superseded by direct filtering on CD8.TEX); 05 composition test and visualization; 06 UCell functional scoring (CD8.TEX); 07 functional scoring patient-level test; 08 functional scoring visualization; 09 CollecTRI pseudobulk TF activity (09b regulon coverage, 09c network coverage); 10 CollecTRI visualization; 11 TCR repertoire clonal analysis.*

- **CD8 purification strategy**: rather than subsetting on the CellTypist-derived `lineage_CD8` flag, CD8 T cells were re-isolated independently on the full TME object of each histology using scGate (`TME_HiRes` model, human), following the same isolation-then-projection strategy used in the official construction of the human CD8 TIL reference atlas (Andreatta et al., 2021, *Nat Commun*; atlas construction pipeline: scGate → UCell cycling/ISG+ removal → STACAS). This avoids circularity with the project's own CellTypist-based annotation and keeps downstream steps methodologically consistent with the reference space used for later ProjecTILs projection.

- **scGate models**: the `TME_HiRes` models (carmonalab/scGate_models) were retrieved with `get_scGateDB()` on 2026-09-17; the scGate package version is recorded in `renv.lock`.
  
- **Performance and memory**: scGate recalculates a dimensionality reduction at each hierarchical gating level by default; supplying the pre-computed Harmony embedding (`reduction = "harmony"`) instead avoids this recomputation and was necessary for tractable runtime on datasets of this size. Memory was a second constraint: on the full LUSC object (457,403 cells), scGate attempted to convert the expression matrix to a dense format, requiring about 96 GB of RAM against 32 GB available. scGate was therefore run in batches of whole patients (about 25,000 cells per batch; 11 batches in LUAD, 19 in LUSC), with patients assigned at random (fixed seed) and response groups spread across batches as evenly as the number of patients allowed. The same batched procedure was applied to both histologies so that no difference between them can be attributed to the way scGate was run, and the batch labels are stored in the object (`scGate_batch`).
  
- **Implication**: the kNN smoothing of scGate scores is computed within each batch, so cells close to the decision threshold can receive a different label than in a whole-object run.
  
- **Validation (LUAD only, the only histology where a whole-object run was feasible)**: the batched run was compared with a whole-object run, for the CD8T population overall, within each response group and per patient. The results are reported in Results by block (Bloc6, "Population definition and cross-validation"). The same comparison could not be performed in LUSC.

- **Cross-check against CellTypist**: two independent identifications of CD8 cells were available, `scGate_multi == "CD8T"` and the CellTypist-derived `lineage_CD8`. The scGate result was carried forward as the working population, to follow the protocol used in building the reference atlas; `lineage_CD8` was kept as a comparison column only. Cells were not filtered on agreement between the two, and discordant cells were left for interpretation after ProjecTILs projection, which uses a CD8-specific reference.

- **Cycling/IFN-high cell removal (UCell)**: applied to the scGate-confirmed CD8T population, using the identical signatures and thresholds from the official CD8 human TIL atlas construction pipeline (cellCycle.G1S, cellCycle.G2M from SignatuR, Tirosh et al., 2016, Science; a manually-defined 4-gene IFN signature, ISG15, IFI6, IFI44L, MX1, which the atlas authors note in their code “works better with UCell than the full signatures from SignatuR package”; script 1-QC_CD8_atlas.Rmd, (https://github.com/carmonalab/CD8_human_TIL_atlas_construction/blob/master/1-QC_CD8_atlas.Rmd)). Thresholds: cellCycle.G1S < 0.1, cellCycle.G2M < 0.1, IFN < 0.25.
  
- **ProjecTILs functional-state projection**: the purified, non-cycling CD8 population was projected onto the official human CD8 TIL reference atlas (Andreatta et al., 2021, *Nat Commun*; atlas construction: `carmonalab/CD8_human_TIL_atlas_construction`), using `Run.ProjecTILs()` with `filter.cells = FALSE` (query already pre-purified, per official guidance for pre-filtered data) and `split.by = "sampleID"` (per-patient batch correction, as recommended).
  
- **Cross-validation and identity correction**: `functional.cluster` (ProjecTILs) was crossed against `annotation_final` (CellTypist) to assess concordance between the two independent methods. Where a substantial group of cells received incompatible labels, identity was settled from raw expression of lineage-defining genes rather than from either label. Mean expression amplitude was the primary criterion, detection rate being treated with caution since ambient RNA can make a gene appear present in cells that do not express it. These cells, labelled as regulatory T cells by CellTypist, showed high CD8A and CD8B with near-absent CD4 and FOXP3, and were counted as CD8 in Bloc6 and Bloc8. Their original labels were not modified, so the Bloc5 composition statistics still count them among regulatory T cells.
  
- **CD8 composition by pathological response**: within the revised CD8 lineage, per-patient proportions of each ProjecTILs functional state were compared between MPR and non-MPR tumors using the Wilcoxon rank-sum test, with Benjamini-Hochberg correction applied across the seven states. Patients with fewer than 50 CD8 cells were excluded from this analysis, since rare states (such as CD8.TPEX, ~1–3% of CD8) have a high probability of a false zero from under-sampling alone below this threshold. This patient-level approach follows Hu et al. (2023, *Genome Medicine*), a comparable neoadjuvant anti-PD-1 NSCLC scRNA-seq study with a similarly small and unbalanced MPR/non-MPR cohort, and supersedes an earlier cell-level chi-squared/`prop.test` approach (see note below). Color coding reused the official ProjecTILs atlas palette; p-value display was standardized using `scales::pvalue()`.
  
- **Note on statistical unit (pseudoreplication)**: cell-level tests (chi-squared, `prop.test`) treat each cell as an independent observation despite within-patient correlation, which can inflate statistical power. This was identified while evaluating `propeller` (Phipson et al., 2022, *Bioinformatics*) on the Bloc5 TME composition data, where a cell-level-significant result did not survive patient-level testing. All composition and functional-scoring comparisons in this project were subsequently re-run at the patient level, including the Bloc5 TME composition analysis itself (see above). That re-analysis is instructive: the cell-level tests had reported an enrichment of macrophages in non-MPR LUAD tumors and an enrichment in MPR LUSC tumors, i.e. an apparent reversal of direction between histologies. At the patient level neither difference is significant, and the medians point in the same direction in both histologies. The reversal was an artefact of treating tens of thousands of cells as independent observations. The pipeline is now patient-level from end to end, with no exception. A formal interaction test (histology × response) was considered but not pursued, given the absence of precedent in comparable published studies, including the original dataset publication (Liu et al., 2025, *Cell*); LUAD/LUSC comparisons are presented descriptively.
  
- **Functional scoring within CD8.TEX (UCell)**: to quantify the intensity of specific functional programs beyond the discrete ProjecTILs states, three short published gene signatures (cytotoxicity: *GZMB, PRF1, FASLG*; exhaustion: *HAVCR2, CTLA4, PDCD1, TIGIT, LAG3*; stemness: *TCF7, SELL, LEF1, IL7R*; Carmona et al., 2020, *Oncoimmunology*, retrieved via `SignatuR$Hs$Programs$Tcell.*`) were scored with UCell within the CD8.TEX population, restricted to the revised CD8 lineage. These signatures were chosen over an alternative source (`signaturesHumanTILs`, ProjecTILs atlas-derived state signatures also available from the Carmona lab) because the latter largely reproduce the discrete states already assigned by ProjecTILs itself, whereas the SignatuR programs provide a continuous, independent readout of functional intensity within those states. Per-patient median scores were compared between MPR and non-MPR using the Wilcoxon rank-sum test with Benjamini-Hochberg correction. Patients contributing fewer than 10 CD8.TEX cells were excluded, following the default filtering used in muscat (Crowell et al., 2020), since a median computed over very few cells is unreliable; this retained 47 LUAD patients (9 MPR, 38 non-MPR) and 89 LUSC patients (25 MPR, 64 non-MPR). Note that this threshold differs from the 50-cell minimum used for the composition analysis because the two apply to different denominators: composition requires enough total CD8 cells to compute a proportion, whereas functional scoring requires enough CD8.TEX cells specifically to compute a reliable median within that state.
  
- **Exclusion of CD8.TPEX from functional analyses**: applying the same 10-cell minimum, only 4 MPR patients would have qualified for CD8.TPEX in LUAD, compared with 13 in LUSC. Analysing this state in LUSC alone would have introduced an asymmetric analytical framework between histologies, so CD8.TPEX was excluded from both functional analyses (signature scoring and transcription factor activity). It remains included in the composition analysis, where the seven states are tested together at the patient level, and in the descriptive tables of cells and patients per state. This focus on the terminally exhausted state is also consistent with its reported clinical relevance: baseline CD8.TEX proportion has been shown to be higher in non-responders than responders prior to immunotherapy and to independently predict checkpoint inhibitor resistance (Fang et al., 2026).
  
- **Transcription factor activity within CD8.TEX (CollecTRI)**: to complement gene-signature scoring with a network-based measure of regulatory activity, transcription factor (TF) activity was inferred using `decoupleR`'s univariate linear model method (`run_ulm`) with the CollecTRI network (Müller-Dott et al., 2023, *Nucleic Acids Research*; ~1,183 human TFs, downloaded directly from static.omnipathdb.org after `decoupleR::get_collectri()` failed due to an unresolved package issue, reformatted to match the expected input, with directionally-ambiguous interactions, 0.4% of the network, excluded). Rather than testing each TF's cell-level activity across patients (which yielded only 5 significant TFs out of 755 tested in an initial analysis, given limited power at n=9-10 MPR patients), a pseudobulk approach was used: raw counts were aggregated per patient within CD8.TEX, a differential expression analysis was run with DESeq2 (MPR vs non-MPR), and `run_ulm` was applied once to the resulting per-gene differential statistic. This approach, recommended in decoupleR's own documentation for this type of comparison, substantially increased sensitivity without introducing repeated per-TF testing across patients. No TFs were pre-selected; all those passing the network's minimum target-gene coverage (`minsize = 5`) were tested.
  
- **Gene filtering harmonised with Bloc7**: `filterByExpr` (Chen et al., 2016) replaces the earlier project-specific filter, which required total counts of at least ten summed across patients and differed from the rule used in the myeloid compartment (see Bloc7). The contrast was also reversed to MPR versus non-MPR, so that a positive score means higher activity in responders in both blocs. The effect on this bloc was substantial, the earlier filter having been the more permissive of the two. With the harmonised filter, 7 factors reach significance in LUAD and none in LUSC.
  
- **CD8.TPEX and TF activity**: consistent with its exclusion from signature scoring, CD8.TPEX was not analysed with CollecTRI.

- **Network coverage check**: because no factor reached significance in LUSC despite a larger differential gene set than in LUAD, the coverage of those genes by the network was examined after the fact, counting the documented regulators of each differentially expressed gene and comparing the distributions between histologies.
  
- **TCR repertoire analysis within CD8.TEX**: to test whether the MPR/non-MPR asymmetry observed in the transcription factor analysis reflects a less-engaged TEX compartment in responders, clonal expansion was assessed using TCR data as a line of evidence independent of gene expression — clone size records past proliferation rather than current transcriptional state. Clonotype assignments were taken directly from the authors' annotation file accompanying GSE243013 (`GSE243013_T_with_TCR_annotation.csv.gz`), joined to the Seurat object by cell barcode; clonotype identifiers are patient-prefixed and non-overlapping across patients (verified: 255,451 clonotypes, none shared). Two complementary metrics were computed per patient within CD8.TEX: (i) the percentage of cells belonging to an expanded clone, using the ≥3-cell threshold defined by the original dataset publication (Liu et al., 2025, *Cell*, STAR Methods), and (ii) Shannon diversity (H = -Σ pᵢ ln pᵢ). Both were calculated explicitly rather than through package defaults, and the pre-computed `expansion` column of the annotation file was deliberately not reused, its underlying threshold being undocumented. The two metrics are non-redundant: the first quantifies the share of the compartment held by proliferating lineages, the second the overall balance of the repertoire.
  
- **Denominator and patient threshold for repertoire metrics**: TCR recovery is incomplete and varies between patients (median 74–85% of TEX cells), and several patients yielded no TCR data at all despite substantial TEX cell numbers (LUSC: 7 patients, including one with 468 TEX cells; LUAD: 2 patients). Using all TEX cells as the denominator would have assigned these patients spurious values of 0% expansion and H = 0. Both metrics were therefore computed over TCR-bearing TEX cells only, and patients were required to contribute at least 10 such cells — the same 10-cell minimum applied in the functional scoring analysis, but applied to the TCR-bearing subset rather than to all TEX cells. This retained 44 LUAD patients (9 MPR, 35 non-MPR) and 81 LUSC patients (24 MPR, 57 non-MPR). Per-patient values were compared between MPR and non-MPR using the Wilcoxon rank-sum test with Benjamini-Hochberg correction across the two metrics, consistent with all other patient-level comparisons in Bloc6.
  
- **Sampling-depth bias and its control**: Shannon diversity is systematically underestimated when few cells are sampled (Mika et al., 2025). Consistent with this, per-patient Shannon diversity increased with the number of cells sampled in the present data (Spearman rho = 0.621 in LUAD, 0.523 in LUSC), so a group sampled less deeply could appear less diverse for technical reasons alone. To rule this out, the number of TCR-bearing CD8.TEX cells per patient was therefore compared between MPR and non-MPR patients before interpreting the diversity results (see TCR repertoire in the results below).

- **Two sensitivity analyses, status fixed before the analysis was run**: metrics recomputed after rarefying patients to a common depth of 30 cells (100 random draws per patient, averaged), which excluded patients with fewer than 30 TCR-bearing cells (37 of 44 retained in LUAD, 70 of 81 in LUSC); and clone size redefined over all T cells rather than within CD8.TEX, since the source threshold does not specify a reference population. The unrarefied within-TEX analysis remains the primary analysis throughout.

### Bloc7_TAM

*Script order (`Scripts/Bloc7_TAM/`): 01 preflight checks; 02 alveolar label check (abandoned, kept for traceability); 03 scGate model inspection; 04 deprecated (superseded by 04b); 04b scGate classification of the myeloid compartment; 04c composition test (Macrophage_pure and Monocyte_pure); 05–10 exploratory functional and subtype scoring, described in `Exploratory_analyses.md`; 11–12 deprecated (perimeter sensitivity, superseded by the Macrophage_pure decision); 13–14 CollecTRI TF activity and overlap; 15–17 DE gene counts, regulon coverage, CollecTRI visualization.*

- **Population definition and classification (scGate)**: rather than treating CellTypist's Confirmed_Myeloid_Bloc5/Unresolved_Bloc5 status as final, cells were classified by scGate, restricted to the myeloid compartment (`lineage_TAM | lineage_TAM_related`) rather than the full object. Unlike Bloc6 (CD8), where scGate had to remain independent of the annotation because ProjecTILs was applied afterward, no reference-based classifier exists for TAM, so scGate itself served as the identity-deciding step, restricted to cells already flagged as myeloid to avoid diluting genuine identity ambiguity with cells that are simply not myeloid. Two models were run, `Macrophage` and `Monocyte` (verified mutually exclusive by construction), classifying cells into four groups: `Macrophage_pure`, `Monocyte_pure`, cells retained by neither model (`scGate_unresolved`), and cells retained by both (`Ambiguous_double_pure`). scGate was adopted as the deciding step because it left substantially fewer cells unresolved than the reference annotation alone in this cohort (13.5% versus 50.2% in LUAD, 23.3% versus 34.6% in LUSC). This is used as an operational criterion, not as evidence that either method is biologically more accurate. `scGate_unresolved` and `Ambiguous_double_pure` cells are reported as a proportion of the compartment but not characterised further and are excluded from downstream analysis: no consensus reference atlas exists for a transitional monocyte-macrophage state in human tumour-infiltrating myeloid cells.

Unlike Bloc6, no cycling or ISG filtering was applied to the macrophage population: that filtering existed to match the ProjecTILs atlas, which excludes such cells by construction, whereas Prolif-TAM and IFN-TAM are recognised TAM states (Ma et al., 2022) and no reference atlas is used here.

| Categories               | LUAD n | LUAD % | LUSC n | LUSC % |
|---------------------------|-------:|-------:|-------:|-------:|
| Confirmed_Myeloid_Bloc5   | 16 607 | 49.5%  | 30 187 | 65.1%  |
| Unresolved_Bloc5          | 16 846 | 50.2%  | 16 043 | 34.6%  |
| Other                     |     77 |  0.2%  |    160 |  0.3%  |
| **Total**                 | **33 530** | **100%** | **46 390** | **100%** |

Cross-tabulating Confirmed_Myeloid_Bloc5/Unresolved_Bloc5 against `scGate_TAM_class` shows that scGate resolves 81.7% (LUAD) and 80.3% (LUSC) of cells left unresolved by CellTypist into `Macrophage_pure` or `Monocyte_pure`, while conversely 16.6% (LUAD) and 28.9% (LUSC) of CellTypist-confirmed myeloid cells are not retained as pure by either scGate model. This divergence reflects the two methods' different decision logics (CellTypist informed by cluster membership, scGate through a hierarchical, per-cell marker gate) rather than an error in either annotation.

- **Composition of the pure myeloid population**: the proportion of `Macrophage_pure` among pure myeloid cells (`Macrophage_pure` and `Monocyte_pure`) was compared between MPR and non-MPR at the patient level with a Wilcoxon rank-sum test, one test per histology (`LUAD_TAM_Proportions`, `LUSC_TAM_Proportions`). The two proportions are complementary, so a single test is sufficient. `scGate_unresolved` and `Ambiguous_double_pure` cells were excluded from the denominator. Patients contributing fewer than 50 pure myeloid cells were excluded, as for CD8 composition, since a proportion computed over a very small denominator is unstable.
  
- **Transcription factor activity (CollecTRI)**: inferred with decoupleR on the CollecTRI network, following the approach applied to CD8.TEX in Bloc6, so that both compartments are analysed in the same way. CollecTRI infers regulatory activity from target gene expression, independently of any cell type signature.
  
- **Two populations analysed separately**: `Macrophage_pure` and `Monocyte_pure`, the two classes of the myeloid compartment retained by scGate (Script04b). Cells classified as `scGate_unresolved` or `Ambiguous_double_pure` are excluded. Separate analysis rather than pooling was chosen so that differences in the balance between monocytes and macrophages across patients could not alter the pseudobulk expression profile of a pooled myeloid compartment. Technical note: the class is stored in the `scGate_TAM_class` column of `Objects/04b_{HIST}_seu_TAM_classified.rds`.
  
- **Design**: pseudobulk per patient (raw counts summed, never averaged normalised data), gene filtering with edgeR's `filterByExpr` (Chen et al., 2016), DESeq2 for the MPR vs non-MPR contrast, then decoupleR `run_ulm` on the per-gene DESeq2 statistic with `minsize = 5`. `filterByExpr` is the recommended strategy for pseudobulk: it keeps genes whose counts per million exceed a threshold in a number of samples set by the design, essentially the size of the smallest group. Using counts per million rather than raw counts avoids favouring patients contributing more cells, which matters here since cell numbers per patient are very unequal. In practice a DGEList is built solely to run the filter, and the retained genes are passed to DESeq2. Network filtered as in Bloc6: complexes excluded, interactions annotated as both stimulatory and inhibitory excluded, mode of regulation encoded as +1 for activation and −1 for repression — note that `is_stimulation` and `is_inhibition` are two separate binary columns, not one signed column.
  
- **Filter harmonisation.** Bloc6 and Bloc7 initially used different gene filters — total counts of at least ten across patients in Bloc6, non-zero in at least a quarter of patients in Bloc7. Neither is wrong, and each comparison between histologies remained valid within its own bloc, but the inconsistency between blocs could not be defended once they are compared. Both were therefore rerun with `filterByExpr`.
  
- **Two thresholds at different levels**: `MIN_CELLS = 10` per patient within each group, so a patient may be retained for one group and excluded from another; and `MIN_PATIENTS = 4` per response group before running the contrast, since with fewer patients the between-patient variability is too imprecisely estimated for a meaningful test. The second threshold is an addition relative to Bloc6, the risk being real with only ten MPR patients in LUAD.
  
- **Reading rule — group without forcing generalisation**: families defined by protein structure (NF-κB, STAT, IRF, C/EBP, nuclear receptors) are established nomenclature and used as such. Functional groupings beyond these require bibliographic support or are dropped. The shared set is not labelled an inflammatory program: CEBPB is reported both as sustaining M1 and as driving M2 expression, IRF1 promotes M1 differentiation while IRF3 inhibits it, and MYC promotes M2 polarisation.
  
- **Bibliographic verification**: performed after the analysis, never before, and restricted first to the myeloid compartment. Only where no myeloid evidence existed were other compartments considered, and this is declared. Pitfall encountered: a search on ZBTB4 returns results on ZBTB14 and ZBTB25, which are different genes, always check the result concerns the right gene.

### Bloc8_Communication

Script order (`Scripts/Bloc8_Communication/`): 01 download of the reference networks; 02 input preparation (per histology); 03 MultiNicheNet analysis and prioritisation (per histology), with 03b geneset assessment under adjusted p-values; 04 filtering of the prioritised interactions (per histology); 05 expression checks (deprecated, not part of the final analysis, kept for traceability only).

- **Scope and rationale**: cell-cell communication was inferred with MultiNicheNet (Browaeys et al., 2023), which builds on NicheNet (Browaeys et al., 2020), restricted to three cell groups: CD8.TEX cells assigned by ProjecTILs (Bloc6) and the scGate-confirmed Macrophage_pure and Monocyte_pure populations (Bloc7). Cells classified as Ambiguous_double_pure or scGate_unresolved were excluded. Macrophage_pure and Monocyte_pure are kept as two separate groups. Nearly half of the transcription factors associated with response in either population were significant in both (31 of 67 in LUAD, 6 of 13 in LUSC), but the communication analysis relies on ligand and receptor expression, and pooling the two populations would make the pseudobulk profile of each patient depend on the balance between monocytes and macrophages in that patient, which may differ between responders and non-responders. Keeping them separate prevents a difference in composition from being read as a difference in expression. Other CD8 states were not included: the analysis is centered on CD8.TEX, and adding groups would increase the number of sender-receiver pairs tested, which is costly with nine MPR patients in LUAD.

- **Reference networks**: the NicheNet v2 prior knowledge networks were taken from two Zenodo records, doi:10.5281/zenodo.7074291 (ligand-target matrix) and doi:10.5281/zenodo.10229222 (annotated ligand-receptor network). They were loaded as in the official MultiNicheNet tutorial: gene aliases converted to official symbols in both objects, then ligand-receptor pairs restricted to ligands present in the matrix, which leaves 5,073 pairs (1,226 ligands, 1,072 receptors) and a matrix of 33,354 target genes by 1,226 ligands.

- **Workflow**: the standard multi-sample workflow was used (pseudobulk differential expression with edgeR through muscat), since every group has at least 9 patients. The sample-agnostic cell-level workflow, designed for datasets with few samples, was not used.

- **Cohort and cell-group sizes**: 

| Histology | Cell group | MPR | non-MPR |
|---|---|---|---|
| LUAD | CD8.TEX | 9 | 38 |
| LUAD | Macrophage_pure | 9 | 34 |
| LUAD | Monocyte_pure | 9 | 38 |
| LUSC | CD8.TEX | 25 | 64 |
| LUSC | Macrophage_pure | 28 | 58 |
| LUSC | Monocyte_pure | 30 | 57 |

Number of patients with at least 10 cells of the cell group (threshold applied per patient and per cell group). Cohort after exclusion of the nine patients treated with chemotherapy alone and of pCR patients: 50 LUAD and 100 LUSC patients.

- **Parameters**: package defaults were kept, since `min_cells = 10` is both the package default and the threshold used across this project. A gene counts as expressed if it has non-zero counts in at least 5 percent of cells of a group in a sample (`fraction_cutoff = 0.05`) and in at least half the samples of the smaller group (`min_sample_prop = 0.50`). Differential expression uses pseudobulk with edgeR, `logFC_threshold = 0.50` and `p_val_threshold = 0.05`, applied to raw p-values (`p_val_adj = FALSE`) in both histologies (see below). The ligand activity criterion of the prioritisation considers upregulated activity only (`ligand_activity_down = FALSE`, the default from version 2 onward), with `top_n_target = 250`.

- **Both contrast directions must be declared explicitly**: MultiNicheNet does not test the reverse comparison automatically, so declaring a single direction would only reveal interactions specific to one group. The syntax is strict and tolerates no whitespace. Note that `make.names()` converts "non-MPR" to "non.MPR", and the contrast labels must match the converted labels exactly or the analysis fails silently; a verification step was built into the preparation script.

- **`p_val_adj = FALSE` in both histologies**: the tutorial uses raw p-values by default and typically recommends adjusted p-values when the smaller group has at least 20 samples, which applies to LUSC (25 MPR patients for CD8.TEX) but not to LUAD (9 MPR patients). Because this is a rule of thumb, the choice was checked with the geneset-to-background ratio that the tutorial prescribes. This ratio divides the number of differentially expressed genes in a cell group by the number of genes considered expressed in it, and should lie between 0.005 and 0.1 (1/200 to 1/10) for the ligand activity step, an enrichment procedure, to be informative: with too few genes no enrichment can be measured, and with too many everything appears enriched. With raw p-values, all six ratios (three cell groups, upregulated and downregulated gene sets) were within range in each histology: 0.024 to 0.069 in LUAD and 0.0085 to 0.031 in LUSC. With adjusted p-values, all six ratios fell below 0.005 in LUAD, with 0 to 44 differentially expressed genes per combination, and three of six fell outside the range in LUSC, where monocytes yielded no gene upregulated in MPR patients and nine upregulated in non-MPR patients. Raw p-values were therefore kept in both histologies, so that the two analyses are handled identically.

- **Prioritisation ranks, it does not test**: the prioritisation score aggregates four criteria weighted equally in the `regular` scenario (up-regulation of the ligand and receptor, cell type specificity of their expression, fraction of samples expressing them, and predicted ligand activity). It carries no significance threshold, which is a departure from the Wilcoxon-based logic used elsewhere in this project. Interactions were assigned to groups with `get_top_n_lr_pairs`, as in the tutorial: each interaction is assigned to the group in which it is prioritised, and interactions whose ligand and receptor are not expressed in that group are discarded. All the interactions passing these steps and the filters described below were exported in full for each histology, with no top-N cut-off. Rank and score are kept in the exported tables but are not used as an argument in the text, since a well-ranked interaction is not necessarily a biologically meaningful one.

- **Structural limitation**: like any communication inference from dissociated single-cell data, MultiNicheNet reads ligand and receptor expression without knowing whether the cells were physically adjacent. Two populations may express a matching pair while occupying distinct tumor regions. The interactions reported are therefore predicted, not established, and spatial transcriptomics would be the validation route.
  
- **Pairs removed before interpretation**: two sets of interactions were removed from the prioritised interactions of each histology. First, interactions whose sender and receiver are the same cell type (referred to as autocrine pairs): they represented 1,867 of 5,326 prioritised interactions in LUAD and 1,633 of 4,666 in LUSC (35% in both), leaving 3,459 and 3,033 interactions respectively. In dissociated cells, a signal sent to the same cell cannot be distinguished from a signal sent to a neighbouring cell of the same type, so these pairs cannot be interpreted as communication from one cell to another. Second, interactions between Macrophage_pure and Monocyte_pure were removed, because the question concerns communication between CD8.TEX cells and the myeloid compartment: only interactions with CD8.TEX as sender or receiver were kept. This left 1,943 interactions in LUAD (655 in MPR, 1,288 in non-MPR) and 1,698 in LUSC (668 in MPR, 1,030 in non-MPR). These two filters are the only modifications made to the selection described in the tutorial.

- **Known axes from the literature**: instead of verifying each interaction individually, ligand–receptor axes described in the literature as immunosuppressive or anti-tumour were identified from published studies, and the exported tables were searched for them. An axis was reported only when the literature supports three elements: the ligand–receptor pair, its production by the sending cell type, and its expression by the receiving cell type. Macrophages and monocytes are not distinguished in this literature validation, a source referring to TAMs, macrophages, monocytes or myeloid cells being taken to cover the myeloid side; the two populations remain analysed separately in the tables. For each axis, the table of axes (given in full in Results by block of Bloc8) lists the effect described by the source, the references, the cancer and model studied, and whether the interaction appears in the prioritisation output of each of the four conditions (LUAD and LUSC, MPR and non-MPR). Presence does not indicate strength or significance. The complete tables are available in the repository `(04_LUAD_filtered_interactions_FULL.csv and 04_LUSC_filtered_interactions_FULL.csv)`

---

## Results by block

### Bloc5_Annotation_TME

**LUAD**

- Every CellTypist_majority subgroup, in each of the 20 clusters, was individually verified against its official marker (Mullan et al., 2023; Domínguez Conde et al., 2021/2022) or, where unavailable, against the trained model's own top marker genes. 51.3% of all cells received a marker-confirmed identity; 48.7% remained explicitly labeled "Unresolved" rather than assigned an unverified identity, and 0.06% were excluded as non-immune (lung epithelial) contamination, a pattern independently corroborated by the CellTypist authors' own documentation of an equivalent artifact.
- `55,841 cells` were confirmed as CD8 lineage (Tem/Trm and Tem/Temra cytotoxic T cells). `26,907 cells` carried a Macrophage/TAM CellTypist label (Macrophages, Intermediate macrophages, Alveolar macrophages, regardless of marker-confirmation status, since transitional states remain part of the biological continuum); `6,623 cells` carried a monocyte label (Classical + Non-classical monocytes), counted separately.

*Patient-level comparison, whole TME (MPR vs non-MPR, Wilcoxon, BH across all 14 displayed categories, patients with ≥50 annotated cells):*

| Population | Median MPR (%) | Median non-MPR (%) | p_value | p_adj |
|---|---|---|---|---|
| CRTAM+ gamma-delta T cells | 1.58 | 0.523 | 0.0971 | 0.404 |
| DC1 | 0.126 | 0.215 | 0.135 | 0.404 |
| Migratory DCs | 0.139 | 0.207 | 0.0873 | 0.404 |
| Regulatory T cells | 8.18 | 11.7 | 0.144 | 0.404 |
| Tcm/Naive helper T cells | 0.754 | 0.307 | 0.0569 | 0.404 |
| CD16+ NK cells | 1.69 | 0.684 | 0.274 | 0.561 |
| Plasma cells | 0.935 | 1.82 | 0.280 | 0.561 |
| Age-associated B cells | 0.246 | 0.187 | 0.452 | 0.585 |
| Classical monocytes | 1.82 | 0.963 | 0.431 | 0.585 |
| DC2 | 1.47 | 1.68 | 0.445 | 0.585 |
| MAIT cells | 0.819 | 0.532 | 0.459 | 0.585 |
| Tem/Trm cytotoxic T cells | 19.2 | 16.5 | 0.676 | 0.789 |
| Tem/Effector helper T cells | 0.638 | 0.681 | 0.858 | 0.924 |
| Macrophages | 0.905 | 0.858 | 0.990 | 0.990 |

No population reaches significance after correction in LUAD. Regulatory T cells show the same direction as in LUSC (higher in non-MPR) but do not survive correction, unlike LUSC where they were the only significant population. Neither Macrophages nor Tem/Trm cytotoxic T cells (the display-level proxy for the CD8 compartment) reach significance; both are retested on their more rigorously defined population in their dedicated blocks.

**LUSC**

- Every `CellTypist_majority` subgroup, in each of the 20 clusters, was individually verified against its official marker (Mullan et al., 2023; Domínguez Conde et al., 2021/2022), independently of its proportion within the parent cluster. 75.2% of all cells received a marker-confirmed identity, 24.8% remained "Unresolved" and 0.03% were labeled as non-immune contamination.

- `102,784 cells` were confirmed as CD8 lineage (Tem/Trm and Tem/Temra cytotoxic T cells; Trm cytotoxic T cells was not confirmed in any LUSC cluster). `22,131 cells` carried a Macrophage/TAM CellTypist label (Macrophages, Intermediate macrophages, Alveolar macrophages — regardless of marker-confirmation status, since transitional states remain part of the biological continuum); `24,259 cells` carried a monocyte label (Classical + Non-classical monocytes), counted separately.

*Patient-level comparison, whole TME (MPR vs non-MPR, Wilcoxon, BH across all 20 categories, patients with ≥50 annotated cells):*

| Population | Median MPR (%) | Median non-MPR (%) | p_value | p_adj |
|---|---|---|---|---|
| Regulatory T cells | 9.81 | 18.0 | 1.28e-5 | **2.56e-4** |
| MAIT cells | 0.794 | 0.478 | 5.78e-3 | 5.78e-2 |
| ILC3 | 0.328 | 0.250 | 2.49e-2 | 1.25e-1 |
| Macrophages | 0.169 | 0.522 | 1.89e-2 | 1.25e-1 |
| DC1 | 0.181 | 0.119 | 5.71e-2 | 2.28e-1 |
| Tem/Effector helper T cells | 0.397 | 0.307 | 1.23e-1 | 4.11e-1 |
| CD16+ NK cells | 1.39 | 1.05 | 2.00e-1 | 5.00e-1 |
| Memory B cells | 16.3 | 10.5 | 1.87e-1 | 5.00e-1 |
| Migratory DCs | 0.114 | 0.156 | 2.41e-1 | 5.36e-1 |
| Plasma cells | 0.925 | 1.61 | 3.27e-1 | 5.95e-1 |
| Tcm/Naive helper T cells | 0.210 | 0.304 | 3.04e-1 | 5.95e-1 |
| Tem/Temra cytotoxic T cells | 0.0386 | 0.0247 | 3.82e-1 | 6.37e-1 |
| CRTAM+ gamma-delta T cells | 0.436 | 0.483 | 4.74e-1 | 6.87e-1 |
| Classical monocytes | 1.40 | 1.38 | 5.49e-1 | 6.87e-1 |
| DC2 | 1.56 | 1.75 | 5.35e-1 | 6.87e-1 |
| Non-classical monocytes | 0.500 | 0.439 | 5.47e-1 | 6.87e-1 |
| Mast cells | 1.58 | 1.39 | 8.68e-1 | 9.64e-1 |
| Tem/Trm cytotoxic T cells | 20.6 | 22.0 | 8.62e-1 | 9.64e-1 |
| Type 1 helper T cells | 0.300 | 0.277 | 9.79e-1 | 9.79e-1 |
| pDC | 0.145 | 0.170 | 9.50e-1 | 9.79e-1 |

### Bloc6_CD8

#### Population definition and cross-validation

**scGate (`TME_HiRes`) classification of the full TME object** — number of cells per scGate label, and share of all cells in each histology.

| scGate label | LUAD (n) | LUAD (%) | LUSC (n) | LUSC (%) |
|---|---|---|---|---|
| CD4T | 69,602 | 26.7 | 130,052 | 28.4 |
| B cell | 58,790 | 22.6 | 81,691 | 17.9 |
| CD8T | 56,347 | 21.6 | 105,302 | 23.0 |
| Macrophage | 19,041 | 7.3 | 12,556 | 2.7 |
| Monocyte | 8,610 | 3.3 | 21,840 | 4.8 |
| NK | 7,244 | 2.8 | 27,860 | 6.1 |
| panDC | 7,511 | 2.9 | 10,809 | 2.4 |
| Plasma cell | 6,635 | 2.5 | 10,179 | 2.2 |
| Mast | 3,936 | 1.5 | 9,338 | 2.0 |
| Neutrophils | 1,499 | 0.6 | 6,693 | 1.5 |
| Epithelial | 4 | <0.1 | 9 | <0.1 |
| Multi (ambiguous) | 2,491 | 1.0 | 4,592 | 1.0 |
| Unclassified (`NA`) | 18,964 | 7.3 | 36,482 | 8.0 |
| **Total** | **260,674** | **100** | **457,403** | **100** |

**Independent scGate isolation agrees moderately with the CellTypist lineage flag**: about seven CD8 cells in ten identified by `lineage_CD8` are also retained by scGate. Concordance is expressed as the share of `lineage_CD8` cells that scGate also classifies as CD8T.

| Histology | scGate CD8T | `lineage_CD8` (CellTypist) | Confirmed by both | scGate only | `lineage_CD8` only | Concordance (% of `lineage_CD8`) |
|---|---|---|---|---|---|---|
| LUAD | 56,347 | 55,841 | 40,312 | 16,035 | 15,529 | 72.2 % |
| LUSC | 105,302 | 102,784 | 72,119 | 33,183 | 30,665 | 70.2 % |

**scGate run in batches of patients (memory constraint) compared with a whole-object run**: because scGate was run in batches of whole patients, the batched result was compared with a whole-object run, which was feasible in LUAD only. The batched run retains 98.6 % of the CD8T cells of the whole-object run (Jaccard index 0.97), with the same agreement in MPR and non-MPR patients and an almost identical ranking of patients by CD8T proportion. The thresholds are practical decision rules set with a margin from the LUAD measurements, not published standards.

| Check | What is tested | Rule | LUAD | LUSC |
|---|---|---|---|---|
| 1a | Patients per response group | Expected counts | MPR 10 / non-MPR 40 (PASS) | MPR 32 / non-MPR 68 (PASS) |
| 1b | Each patient in exactly one batch | Maximum = 1 | 1 (PASS) | 1 (PASS) |
| 1c | Every cell has a batch and a unique name | 0 missing, 0 duplicated | 0 ; 0 (PASS) | 0 ; 0 (PASS) |
| 2a | Whole-run CD8T cells still CD8T in the batched run | ≥ 0.95 | 0.986 (PASS) | not feasible |
| 2b | Jaccard index of CD8T | ≥ 0.95 | 0.970 (PASS) | not feasible |
| 3a | Jaccard index of CD8T within each response group | Difference between groups ≤ 0.02 | MPR 0.976 / non-MPR 0.969 (PASS) | not feasible |
| 3b | Spearman correlation of per-patient CD8T proportion | ≥ 0.99 | 0.9988 (PASS) | not feasible |
| 3c | Largest per-patient change in CD8T count | Within ±10 % | −8.7 % to +2.1 % (PASS) | not feasible |

Checks 2a to 3c require a whole-object run, which does not fit in 32 GB of RAM for LUSC. One LUAD batch (batch 11) contains no MPR cell (10 MPR patients for 11 batches); in LUSC, all 19 batches contain both response groups.

- **Cycling and IFN-high cells were removed to match the ProjecTILs atlas**, which excludes them by construction. **The removal was much heavier in non-MPR than in MPR patients in both histologies** (LUAD: 16.7 % against 4.8 % of scGate CD8T cells; LUSC: 14.7 % against 5.3 %), although the same thresholds were applied to every cell. The IFN criterion dominated (~70 % of cumulative flags in LUAD, ~71 % in LUSC), plausibly reflecting an active immune response under neoadjuvant anti-PD-1, though this is not interpreted biologically here.

| | LUAD | LUSC |
| --- | --- | --- |
| Cells in the object (all scGate labels) | 260,674 | 457,403 |
| scGate CD8T cells before filtering | 56,347 (21.6 %) | 105,302 (23.0 %) |
| Cells after filtering | 47,996 | 93,584 |
| Cells removed (% of scGate CD8T) | 8,351 (14.8 %) | 11,718 (11.1 %) |
| MPR: removed / before (%) | 437 / 9,015 (**4.8 %**) | 2,078 / 39,533 (**5.3 %**) |
| non-MPR: removed / before (%) | 7,914 / 47,332 (**16.7 %**) | 9,640 / 65,769 (**14.7 %**) |
| Patients retained (MPR / non-MPR) | 10 / 40 | 32 / 68 |
| Flagged: cellCycle.G1S ≥ 0.1 | 1,495 | 1,987 |
| Flagged: cellCycle.G2M ≥ 0.1 | 1,355 | 1,933 |
| Flagged: IFN ≥ 0.25 | 6,753 | 9,554 |

Flags are not exclusive: a cell can exceed several thresholds, so their sum (9,603 in LUAD, 13,474 in LUSC) is larger than the number of cells removed.

**ProjecTILs (`filter.cells = FALSE`, projection per patient) assigned each purified CD8 cell to one of the atlas's seven functional states.** In LUAD, CD8.EM and CD8.TEX dominate at similar levels (32.1 % and 31.0 %); in LUSC, CD8.CM is the largest state (33.8 %) and CD8.TEX is lower (23.1 %). CD8.TPEX remains rare in both histologies (3.0 % and 1.4 %) and is absent in a number of patients (7/50 in LUAD, 15/100 in LUSC), consistent with it not being testable in the functional analyses (see below).

**ProjecTILs functional states of the purified CD8 population**

| State | LUAD cells | LUAD % | LUSC cells | LUSC % |
|---|---|---|---|---|
| CD8.EM | 15,423 | 32.1 | 27,063 | 28.9 |
| CD8.TEX | 14,865 | 31.0 | 21,571 | 23.1 |
| CD8.CM | 12,542 | 26.1 | 31,602 | 33.8 |
| CD8.TEMRA | 1,611 | 3.4 | 6,556 | 7.0 |
| CD8.TPEX | 1,429 | 3.0 | 1,341 | 1.4 |
| CD8.MAIT | 1,133 | 2.4 | 2,560 | 2.7 |
| CD8.NaiveLike | 993 | 2.1 | 2,891 | 3.1 |
| **Total** | **47,996** | **100** | **93,584** | **100** |

Patients per state (of 10 MPR / 40 non-MPR in LUAD, 32 MPR / 68 non-MPR in LUSC):

| State | LUAD MPR | LUAD non-MPR | LUSC MPR | LUSC non-MPR |
|---|---|---|---|---|
| CD8.CM | 10 | 39 | 32 | 67 |
| CD8.EM | 10 | 39 | 32 | 68 |
| CD8.MAIT | 10 | 37 | 31 | 59 |
| CD8.NaiveLike | 10 | 35 | 29 | 60 |
| CD8.TEMRA | 9 | 36 | 32 | 57 |
| CD8.TEX | 10 | 39 | 28 | 67 |
| CD8.TPEX | 10 | 33 | 26 | 59 |

**Per-patient testability of ProjecTILs states**: share of patients with at least one cell in each state, and, among those, share with fewer than 10 cells — the threshold used for UCell functional scoring within a state.

| State | LUAD ≥1 cell (%) | LUAD, of those <10 cells (%) | LUSC ≥1 cell (%) | LUSC, of those <10 cells (%) |
|---|---|---|---|---|
| CD8.CM | 98 | 2.0 | 99 | 2.0 |
| CD8.EM | 98 | 2.0 | 100 | 3.0 |
| CD8.MAIT | 94 | 31.9 | 90 | 38.9 |
| CD8.NaiveLike | 90 | 42.2 | 89 | 40.4 |
| CD8.TEMRA | 90 | 28.9 | 89 | 27.0 |
| CD8.TEX | 98 | 4.1 | 95 | 6.3 |
| CD8.TPEX | 86 | 39.5 | 85 | 54.1 |

- **CD8.TPEX was not tested in the functional analyses of either histology.** With the 10-cell minimum, only 4 MPR patients in LUAD retained this state (13 in LUSC), too few for a meaningful MPR versus non-MPR comparison in LUAD. Testing it in LUSC alone would have created an asymmetric framework between histologies, so it was excluded from both. CD8.TEX, by contrast, met the threshold in most patients (9 of 10 MPR in LUAD, 25 of 32 in LUSC) and is analysed throughout.

#### Composition of functional states

- **The baseline profile differs between histologies.** LUAD is dominated by CD8.EM (32.1 %), CD8.TEX (31.0 %) and CD8.CM (26.1 %); LUSC by CD8.CM (33.8 %), CD8.EM (28.9 %) and CD8.TEX (23.1 %) — a LUSC profile shifted toward less-differentiated, memory-like states and away from terminal exhaustion.

- **Composition distinguishes response groups in LUSC only.** At the patient level (Wilcoxon, BH across the seven states; patients with fewer than 50 CD8 cells excluded — 48/50 retained in LUAD, 96/100 in LUSC), CD8.TEX was higher in non-MPR LUSC tumors (median 22.6 versus 9.0 %; adjusted P = 0.005), with CD8.CM (35.8 versus 24.9 %; adjusted P = 0.007) and CD8.MAIT (3.1 versus 1.2 %; adjusted P = 0.007) higher in MPR. In LUAD, no state reached significance (all adjusted P ≥ 0.89).

LUAD

| State | Median MPR (%) | Median non-MPR (%) | p_adj |
|---|---|---|---|
| CD8.CM | 29.5 | 26.7 | 0.889 |
| CD8.EM | 24.4 | 29.0 | 0.889 |
| CD8.MAIT | 2.0 | 2.2 | 0.889 |
| CD8.NaiveLike | 2.2 | 1.4 | 0.889 |
| CD8.TEMRA | 4.7 | 1.9 | 0.889 |
| CD8.TEX | 21.9 | 17.6 | 0.889 |
| CD8.TPEX | 2.0 | 1.5 | 0.889 |

LUSC

| State | Median MPR (%) | Median non-MPR (%) | p_adj |
|---|---|---|---|
| CD8.CM | 35.8 | 24.9 | **0.007** |
| CD8.EM | 22.1 | 26.8 | 0.215 |
| CD8.MAIT | 3.1 | 1.2 | **0.007** |
| CD8.NaiveLike | 1.6 | 1.5 | 0.808 |
| CD8.TEMRA | 2.6 | 2.0 | 0.215 |
| CD8.TEX | 9.0 | 22.6 | **0.005** |
| CD8.TPEX | 0.9 | 1.0 | 0.543 |

#### Functional scoring within CD8.TEX

- **The same pattern holds in both histologies, and it is not a simple gradient.** Patients with fewer than 10 CD8.TEX cells were excluded (LUAD: 47/48 retained, 9 MPR / 38 non-MPR; LUSC: 89/96 retained, 25 MPR / 64 non-MPR). CD8.TEX cells in non-MPR tumors were both more cytotoxic and more exhausted, in both histologies. The two are not mutually exclusive, exhausted CD8 cells being known to retain partial cytotoxic potential.

| Signature | LUAD MPR | LUAD non-MPR | LUAD p_adj | LUSC MPR | LUSC non-MPR | LUSC p_adj |
|---|---|---|---|---|---|---|
| Cytotoxicity | 0.284 | 0.393 | **0.008** | 0.310 | 0.400 | **0.001** |
| Exhaustion | 0.184 | 0.314 | **0.007** | 0.225 | 0.304 | **0.001** |
| Stemness | 0 | 0 | 0.665 | 0 | 0 | 0.522 |

- **Stemness did not differ between response groups** (adjusted P = 0.665 in LUAD, 0.522 in LUSC), with a median of zero in all four groups, consistent with cells already committed to terminal exhaustion.

#### Composition versus functional intensity

**The cytotoxicity and exhaustion differences reached significance in both histologies**, whereas the composition analysis reached significance in LUSC only. This contrast should not be over-read: the two analyses draw on different patient subsets (≥50 CD8 cells, 48/50 in LUAD and 96/100 in LUSC, versus ≥10 CD8.TEX cells, 47/48 in LUAD and 89/96 in LUSC), correct across families of different sizes (seven states versus three programs), and were never compared formally. What can be said is that in LUSC the two converge — CD8.TEX are both more abundant and more exhausted in non-MPR tumors — while in LUAD only the functional signal is detected.

#### Transcription factor activity in CD8.TEX (CollecTRI, pseudobulk)

**Regulon activity was inferred from a pseudobulk differential signature (DESeq2, MPR vs non-MPR) run through CollecTRI (`run_ulm`, minsize = 5), on the same patients as the functional scoring above** (LUAD: 47/50 with ≥10 CD8.TEX cells; LUSC: 89/100). A positive score means higher activity in MPR tumors.

| | LUAD | LUSC |
| --- | --- | --- |
| Genes after `filterByExpr` | 7,960 / 34,115 | 7,414 / 34,931 |
| TFs tested (≥5 targets) | 477 | 459 |
| Significant (adjusted P < 0.05) | 7 | 0 |

**In LUAD, STAT1 and RELA had the strongest negative scores** (higher activity in non-MPR tumors), alongside NFKB1, CEBPZ, HSF1, STAT3 and E2F4. **In LUSC, no TF reached significance after correction** (lowest adjusted P = 0.331).

- **Among the factors reaching significance, all lie on the non-responder side in LUAD.** Of 477 factors tested in LUAD, 7 differed significantly, all more active in non-responders: STAT1, RELA, NFKB1, CEBPZ, HSF1, STAT3 and E2F4. None reached significance on the responder side. In LUSC, none of the 459 factors tested reached significance.

**A note on interpretation.** Significance here reflects a detectable difference between MPR and non-MPR, not the presence or absence of TF activity in CD8.TEX cells: a non-significant TF may still be active in this compartment, and the same caveat applies to LUAD's non-significant factors as to LUSC's absence of any significant result.

**A check for whether LUSC's absence of significant TF reflects a lack of differential signal or a small number of regulons shows neither.** Applying the same DESeq2 pseudobulk model at the gene level (before CollecTRI aggregation), 258 genes differed between response groups in LUSC at adjusted P < 0.05 (median |Wald statistic| 0.94), and the number of regulons tested was comparable to that in LUAD (459 in LUSC against 477, with a median size of 15 targets in both). The absence of significant TF in LUSC is therefore not explained by a lack of differential signal or by the number of regulons tested. One possibility is that the LUSC signal is more diffuse across genes and regulons, rather than concentrated on a few dominant TF as in LUAD (where STAT1 and RELA carry the strongest scores), but this is not tested directly here.

**A likely explanation lies in the network coverage of the specific differential genes rather than in regulon-level coverage** Although CollecTRI's overall regulon coverage is comparable between histologies (Script 09b), the genes carrying the differential signal are covered very unevenly.

| | LUAD | LUSC |
| --- | --- | --- |
| Differential genes (adjusted P < 0.05) | 65 | 258 |
| With no regulator in CollecTRI | 24 (36.9 %) | 146 (56.6 %) |
| Median regulators per gene | 2 | 0 |

**TIGIT illustrates the network's blind spots beyond any single histology**: differentially expressed in both LUAD (P = 0.033) and LUSC (P = 0.035), it has no documented regulator in CollecTRI in either. Well-covered genes such as GZMB, SLC2A3, JUNB and BTG2 (14–20 regulators each) are comparatively rare among the differential set.

**Ranked by P value regardless of significance, the leading factors overlap little between histologies, although no LUSC factor reaches significance** (the seven leading factors on each side are given in the tables of Results/Bloc6_CD8/Tables/). On the non-MPR side, LUAD is led by STAT1, RELA, NFKB1, CEBPZ, HSF1, STAT3 and E2F4, which are its seven significant factors; the leading LUSC factors are IRF3, E2F4, STAT1, IRF2, PPARA, IRF1 and IRF9, of which only STAT1 and E2F4 are shared with LUAD. On the MPR side, where no factor reached significance in either histology, the leading factors are ZBED1, ZBTB4, HIVEP2, KAT2B, MSX1, GFI1 and KLF10 in LUAD, and NR2F1, MSX1, RORC, KAT7, HOXA7, TGFB1I1 and ELF3 in LUSC; only MSX1 is shared.

#### TCR repertoire

- **Neither clonal expansion nor diversity differs**, in either histology, under any of the analytical variants tested.

| Histology | Metric | Median MPR | Median non-MPR | BH-adjusted P |
|---|---|---|---|---|
| LUAD | % cells in expanded clones | 66.5 | 65.9 | 0.474 |
| LUAD | Shannon diversity | 3.53 | 3.33 | 0.474 |
| LUSC | % cells in expanded clones | 74.7 | 67.8 | 0.515 |
| LUSC | Shannon diversity | 2.89 | 2.90 | 0.481 |

Patients contributing fewer than 10 TCR-bearing CD8.TEX cells were excluded (LUAD: 44/48 retained, 9 MPR / 35 non-MPR; LUSC: 81/96, 24 MPR / 57 non-MPR). Sensitivity analyses were likewise non-significant: Gini-Simpson (LUAD P = 0.098, LUSC P = 0.22), rarefaction to 30 cells (LUAD adjusted P = 0.148 for both metrics, LUSC 0.394 for both), Gini-Simpson rarefied (LUAD P = 0.117, LUSC P = 0.424), and clone size defined over all T cells rather than CD8.TEX alone (LUAD P = 0.352, LUSC P = 0.653).

Sampling depth, the number of TCR-bearing CD8.TEX cells per patient, did not differ between response groups in either histology:

| Histology | Patients (MPR / non-MPR) | Median TCR-bearing CD8.TEX cells, MPR | Median TCR-bearing CD8.TEX cells, non-MPR | Wilcoxon P |
|---|---|---|---|---|
| LUAD | 9 / 35 | 187 | 168 | 0.67 |
| LUSC | 24 / 57 | 146 | 129 | 0.80 |

Source: `Results/Bloc6_CD8/Tables/11_LUAD_TCR_sampling_depth.csv` and `Results/Bloc6_CD8/Tables/11_LUSC_TCR_sampling_depth.csv`.

- **The median differences should not be read as trends.** For % expanded clones, MPR is nominally higher than non-MPR in both histologies, the opposite of the prediction stated before the analysis, but the difference is negligible in LUAD (+0.6 points) and neither reaches significance. For Shannon diversity, LUAD shows a small difference in the predicted direction (MPR higher) while LUSC shows essentially none (a difference of 0.01). No pattern is consistent enough, in either metric or either histology, to read as a real effect rather than noise.

- **This tests a hypothesis and does not support it.** The idea that the TEX compartment is less engaged in responders was advanced to account for the asymmetry in transcription factor activity. The clonal repertoire does not support it, and that asymmetry therefore remains unexplained.

#### Summary

Across the three planned analyses of the CD8 lineage, functional scoring separated response groups in both histologies, whereas composition reached significance in LUSC only and transcription factor activity in LUAD only.

A fourth analysis, TCR clonality, was added to test whether the asymmetry towards non-responders observed in the transcription factor analysis could reflect reduced clonal engagement of the TEX compartment in responders. It did not: clonal expansion and diversity were indistinguishable between response groups in both histologies, and the asymmetry remains unexplained.

### Bloc7_TAM

Before TAM composition and functional scoring, myeloid coverage was checked per patient. Of the 150 analysed patients (50 LUAD, 100 LUSC), 10 patients (2 LUAD, 8 LUSC) had no myeloid cells (neither TAM nor TAM-related lineages) captured by CellTypist annotation; their immune compartment was entirely lymphoid (T, NK, B, plasma cells). This reflects upstream cell capture, not an artifact of the TAM gating strategy.

| Histology | Analysed patients | Patients with ≥1 myeloid cell | Patients with no myeloid cells | Response of excluded patients |
|---|---|---|---|---|
| LUAD | 50 | 48 | 2 | 2 non-MPR |
| LUSC | 100 | 92 | 8 | 7 non-MPR, 1 MPR |

Within `lineage_TAM | lineage_TAM_related` (Macrophages, Classical monocytes and Non-classical monocytes confirmed by CellTypist, plus all unresolved myeloid labels, including Intermediate macrophages), the myeloid compartment totals 33,530 cells in LUAD and 46,390 cells in LUSC.

| Category                  | LUAD n | LUAD % | LUSC n | LUSC % |
|----------------------------|-------:|-------:|-------:|-------:|
| Confirmed_Myeloid_Bloc5    | 16 607 | 49.5%  | 30 187 | 65.1%  |
| Unresolved_Bloc5           | 16 846 | 50.2%  | 16 043 | 34.6%  |
| Other                      |     77 |  0.2%  |    160 |  0.3%  |
| **Total**                  | **33 530** | **100%** | **46 390** | **100%** |

Within the resulting myeloid compartment (`lineage_TAM | lineage_TAM_related`), cells were cross-tabulated by CellTypist confirmation status (`Confirmed_Myeloid_Bloc5` / `Unresolved_Bloc5` / `Other`) against `scGate_TAM_class` (`Macrophage_pure` / `Monocyte_pure` / `scGate_unresolved` / `Ambiguous_double_pure`).

**LUAD** (33,530 cells)

| CellTypist status (Bloc5) | scGate: Macrophage_pure | scGate: Monocyte_pure | scGate: Ambiguous_double_pure | scGate: unresolved | Total |
|---|---:|---:|---:|---:|---:|
| Confirmed_Myeloid_Bloc5 | 9 980 | 3 778 | 100 | 2 749 | 16 607 |
| Unresolved_Bloc5 | 8 992 | 4 767 | 1 356 | 1 731 | 16 846 |
| Other | 29 | 4 | 3 | 41 | 77 |
| **Total** | **19 001** | **8 549** | **1 459** | **4 521** | **33 530** |

**LUSC** (46,390 cells)

| CellTypist status (Bloc5) | scGate: Macrophage_pure | scGate: Monocyte_pure | scGate: Ambiguous_double_pure | scGate: unresolved | Total |
|---|---:|---:|---:|---:|---:|
| Confirmed_Myeloid_Bloc5 | 5 354 | 14 793 | 1 309 | 8 731 | 30 187 |
| Unresolved_Bloc5 | 5 659 | 7 218 | 1 186 | 1 980 | 16 043 |
| Other | 22 | 37 | 18 | 83 | 160 |
| **Total** | **11 035** | **22 048** | **2 513** | **10 794** | **46 390** |

`Confirmed_Myeloid_Bloc5` aggregates three CellTypist labels (Macrophages, Classical monocytes, Non-classical monocytes) rather than macrophages alone. The scGate columns reflect an independent, marker-based classification of the same cells: a cell counted as `Confirmed_Myeloid_Bloc5` is not necessarily classified `Macrophage_pure` by scGate, and vice versa, since the two methods do not use the same criteria.

Each row shows how cells sharing a given CellTypist status (`Confirmed_Myeloid_Bloc5`, `Unresolved_Bloc5`, `Other`) distribute across the four scGate classes: for example, of the 16,607 `Confirmed_Myeloid_Bloc5` cells in LUAD, 9,980 are classified `Macrophage_pure` by scGate, 3,778 `Monocyte_pure`, 100 `Ambiguous_double_pure`, and 2,749 remain `scGate_unresolved`. `Confirmed_Myeloid_Bloc5` itself aggregates three CellTypist labels (Macrophages, Classical monocytes, Non-classical monocytes); scGate classification is independent of which of the three a cell was assigned.

- **Composition of the pure myeloid population does not differ between response groups**. Patients with at least 50 pure myeloid cells were retained: 44 in LUAD (8 MPR, 36 non-MPR) and 83 in LUSC (28 MPR, 55 non-MPR). The median proportion of Macrophage_pure was 37.4% in MPR and 40.8% in non-MPR tumors in LUAD (Wilcoxon P = 0.777), and 22.3% and 33.1% in LUSC (P = 0.232). Monocyte_pure is the complement of Macrophage_pure, so it carries the same P value.

- **Exploratory functional and subtype scoring of the macrophage compartment** (Cheng et al. signatures, Wang et al. markers, Hallmark metabolic programs) is described and reported separately in `Exploratory_analyses.md`; these analyses are not part of the preprint.

#### Transcription factor activity in the myeloid compartment (CollecTRI, pseudobulk)

**Regulon activity was inferred separately in Macrophage_pure and Monocyte_pure from a pseudobulk differential signature (DESeq2, MPR vs non-MPR) run through CollecTRI (`run_ulm`, minsize = 5).** A positive score means higher activity in MPR tumors.

**Patients retained (pseudobulk, CollecTRI)**

| | LUAD | LUSC |
|---|---|---|
| Macrophage_pure | 9 MPR vs 34 non-MPR | 28 MPR vs 58 non-MPR |
| Monocyte_pure | 9 MPR vs 38 non-MPR | 30 MPR vs 57 non-MPR |

**Transcription factors tested and significant, by histology and population**

| | LUAD Macrophage_pure | LUAD Monocyte_pure | LUSC Macrophage_pure | LUSC Monocyte_pure |
|---|---|---|---|---|
| TFs evaluated | 552 | 512 | 475 | 509 |
| Significant TFs | 64 | 34 | 8 | 11 |
| Higher in MPR | 2 | 1 | 1 | 3 |
| Higher in non-MPR | 62 | 33 | 7 | 8 |

The factors reaching significance are given in the tables of `Results/Bloc7_TAM/Tables/`. In LUAD, 31 factors are significant in both populations, of which 30 are more active in non-MPR tumors and one (ZBTB4) in MPR tumors, with RELA, NFKB1 and STAT1 among the strongest scores in Macrophage_pure; 33 further factors are specific to Macrophage_pure and 3 to Monocyte_pure (PPARA, MYC, SMAD3). In LUSC, 6 factors are significant in both populations (STAT1, RELA, NFKB1, IRF1, REL, TP73), 2 are specific to Macrophage_pure (CIITA, HSF4) and 5 to Monocyte_pure (HIF1A, ZBTB4, IRF9, TRERF1, STAT2).

- **Among the factors reaching significance, most are more active in non-responders, in both histologies and both populations**: 62 of 64 in Macrophage_pure and 33 of 34 in Monocyte_pure in LUAD, 7 of 8 and 8 of 11 in LUSC.

**Differential gene signal by group and histology**

Gene-level counts from the same pseudobulk model, before CollecTRI aggregation, are given for reference. Patients are those retained for each population.

| Histology | Group | Patients | Genes tested (filterByExpr) | padj < 0.05 | padj < 0.01 | median\|stat\| |
|---|---|---|---|---|---|---|
| LUAD | Macrophage_pure | 43 | 9,483 | 14 | 1 | 0.813 |
| LUAD | Monocyte_pure | 47 | 7,961 | 221 | 40 | 1.01 |
| LUSC | Macrophage_pure | 86 | 5,744 | 209 | 35 | 0.991 |
| LUSC | Monocyte_pure | 87 | 7,929 | 27 | 1 | 0.895 |

### Bloc8_Communication

- **Input of the analysis (patients, cells and genes)**: The first table gives, for each histology and cell group, the number of patients passing the 10-cell threshold and the number of cells in the input object. The second gives the number of genes before and after expression filtering.

**Input of the MultiNicheNet analysis**

| Histology | Cell group | Patients MPR | Patients non-MPR | Cells MPR | Cells non-MPR |
|---|---|---|---|---|---|
| LUAD | CD8.TEX | 9 | 38 | 2,232 | 12,633 |
| LUAD | Macrophage_pure | 9 | 34 | 1,035 | 17,965 |
| LUAD | Monocyte_pure | 9 | 38 | 1,213 | 7,336 |
| LUSC | CD8.TEX | 25 | 64 | 5,950 | 15,617 |
| LUSC | Macrophage_pure | 28 | 58 | 4,795 | 6,240 |
| LUSC | Monocyte_pure | 30 | 57 | 10,764 | 11,284 |

**Genes before and after expression filtering**

| Histology | Genes in the input | Genes retained after filtering |
|---|---|---|
| LUAD | 34,115 | 11,836 |
| LUSC | 34,931 | 11,184 |

*Legend:* Patients are those contributing at least 10 cells of the cell group (min_cells = 10, applied per patient and per cell group). Cell counts are those of the input object, before this patient filter. A cell present in both the CD8 object (Bloc 6) and the TAM object (Bloc 7) cannot carry two labels and was excluded: 4 cells in LUAD, of which one was Macrophage_pure (the other three belonging to none of the three cell groups), and 12 cells in LUSC, of which four were CD8.TEX (the other eight belonging to none of the three cell groups). The cell counts above are given after this exclusion. Retained genes are those considered expressed in at least one of the three cell groups (fraction_cutoff = 0.05, min_sample_prop = 0.50).

- **Interactions in the prioritisation output**
After removal of same-cell-type pairs and of macrophage–monocyte pairs (see Methods), the prioritisation output comprised 1,943 interactions in LUAD (655 assigned to MPR and 1,288 to non-MPR) and 1,698 in LUSC (668 and 1,030), corresponding to 820 and 732 distinct ligand–receptor pairs. The number of interactions assigned to each response group is not interpreted as a biological result, the two groups differing in size.

- **Known axes found in the output**
Nine axes described in the literature were found in the output; they are listed with their sources, cancers and models in the table (below). Four are described as immunosuppressive (CD274–PDCD1, CD80/CD86–CTLA4, NECTIN2/PVR–TIGIT and SPP1–CD44), two as anti-tumour (CD80/CD86–CD28 and CXCL9/10/11–CXCR3), and three as context-dependent (CXCL16–CXCR6, CD200–CD200R1 and IFNG–IFNGR1/2). In non-responders, all nine were found in LUAD and eight in LUSC, CXCL16–CXCR6 being the only one absent there. In responders, only CXCL16–CXCR6 was found in LUAD, and CXCL16–CXCR6 and SPP1–CD44 in LUSC; the other seven axes were found in non-responders only, in both histologies. Seven of the nine run from myeloid cells towards CD8.TEX; CD200–CD200R1 and IFNG–IFNGR1/2 run from CD8.TEX towards myeloid cells. Presence in the output means that the interaction was prioritised for that group, not that it is more active there, and differences in the number of axes between groups are descriptive.

**Ligand–receptor axes described in the literature and found in the CD8.TEX–myeloid cell prioritisation output**

| Axis | Direction | Reported effect | Source | Cancer and model | LUAD MPR | LUAD non-MPR | LUSC MPR | LUSC non-MPR |
|---|---|---|---|---|---|---|---|---|
| CXCL16–CXCR6 | Myeloid cells to CD8.TEX | Context-dependent: CXCL16 on myeloid cells recruits and retains CXCR6⁺ CD8 T cells, while contact with CXCL16⁺ myeloid cells inside the tumour has also been linked to T cell dysfunction | Li et al. 2026; Chia et al. 2024 | Review, solid tumours; human and mouse glioblastoma | Present | Present | Present | Absent |
| IFNG–IFNGR1/2 | CD8.TEX to myeloid cells | Context-dependent: IFNγ, of which CD8 T cells are key producers, acts on IFNGR1, expressed by all nucleated cells; initial exposure is described as antitumorigenic, including recruitment of T cells through CXCL9/10/11 and priming of macrophages towards an antitumorigenic phenotype, whereas prolonged exposure promotes protumorigenic effects, including immunosuppression through PD-L1 and IDO1 | Gocher et al. 2022 | Review, several tumour types | Absent | Present | Absent | Present |
| CD80/CD86–CD28 | Myeloid cells to CD8.TEX | Anti-tumour: co-stimulation that CD8 T cells need to be fully activated; without it, they become anergic | Sato et al. 2022; Raskov et al. 2021 | Human NSCLC; review | Absent | Present | Absent | Present |
| CD80/CD86–CTLA4 | Myeloid cells to CD8.TEX | Immunosuppressive: CTLA-4 competes with CD28, dampens TCR signalling and promotes exhaustion; the authors propose that CTLA-4⁺ T cells interact with CD80 on macrophages | Sato et al. 2022; Zhong et al. 2025 | Human NSCLC; review | Absent | Present | Absent | Present |
| NECTIN2/PVR–TIGIT | Myeloid cells to CD8.TEX | Immunosuppressive: TIGIT delivers inhibitory signals, reducing activation and cytokine production and promoting CD8 T cell exhaustion | Zhu et al. 2025; Zhong et al. 2025; Liu et al. 2024 | Human NSCLC; general review; NSCLC review | Absent | Present | Absent | Present |
| CD274–PDCD1 | Myeloid cells to CD8.TEX | Immunosuppressive: PD-L1 bound to PD-1 on exhausted CD8 T cells inhibits TCR signalling, cytotoxicity, cytokine production and proliferation | Shinchi et al. 2022; Zhong et al. 2025; Wei et al. 2024; Sato et al. 2026 | Human LUAD; review; human pan-cancer; human NSCLC | Absent | Present | Absent | Present |
| SPP1–CD44 | Myeloid cells to CD8.TEX | Immunosuppressive: SPP1 secreted by TAMs engages CD44 on T cells; described as a major immunosuppressive pathway in NSCLC, and SPP1⁺ TAMs are strongly linked to resistance to immune checkpoint inhibitors | Zhu et al. 2025; Saldana et al. 2026 | Human NSCLC, scRNA-seq, spatial transcriptomics and in vitro co-culture; review | Absent | Present | Present | Present |
| CXCL9/10/11–CXCR3 | Myeloid cells to CD8.TEX | Anti-tumour: CXCL9/10/11 produced by TAMs are described as attracting CXCR3⁺ CD8 T cells; CXCL9⁺/CXCL10⁺ TAMs are linked to better responses to immune checkpoint inhibitors | Saldana et al. 2026; House et al. 2020; Wei et al. 2024 | Review, NSCLC; mouse models and patients under dual PD-1/CTLA-4 blockade; human pan-cancer | Absent | Present | Absent | Present |
| CD200–CD200R1 | CD8.TEX to myeloid cells | Context-dependent: CD200, carried by activated T cells, engages CD200R1 on myeloid cells and dampens their activation; the net effect on T cells is indirect and described as dependent on context; in NSCLC, CD200R1 is upregulated on TAMs compared with macrophages from normal lung | Liu et al. 2020; Sahu et al. 2026; Tøndell et al. 2021 | Review; review, solid tumours; human NSCLC | Absent | Present | Absent | Present |

*Present: the interaction appears in the MultiNicheNet prioritisation output for that condition; Absent: it does not. Presence does not indicate strength or significance. Items in "Cancer and model" follow the order of the references in "Source". The full interaction tables (`04_LUAD_filtered_interactions_FULL.csv` and `04_LUSC_filtered_interactions_FULL.csv`) are available in the repository.*

---

## Author

Myriam Yasmina Soumahoro, MSc Biology
Independent researcher

ORCID: [0009-0009-8549-5340](https://orcid.org/0009-0009-8549-5340)

Contact: m.yasminasoumahoro@gmail.com