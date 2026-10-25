## Bloc7_TAM: functional and subtype scoring frameworks

### Methodology

These three frameworks (the Cheng et al. functional signatures, the Wang et al. subtype markers, and the Hallmark metabolic programs) are applied here through general-purpose per-cell scoring tools (UCell), independently of the clustering pipelines the source studies used to derive and characterize them.

#### Scope

All three frameworks are applied to scGate-confirmed macrophages (Macrophage_pure), as classified in Script 04b. Classical monocytes (`Monocyte_pure`) were not scored: they have not engaged macrophage differentiation, so the Cheng et al. (2021) signatures would measure states they do not possess.

#### Functional scoring (Cheng et al.)

Four published signatures defining M1, M2, angiogenesis and phagocytosis phenotypes in tumor-infiltrating myeloid cells (Cheng et al., 2021, *Cell*, Table S5) were scored with UCell. Discrete subtype assignment was deliberately not attempted: Cheng et al. show that markers of pro-angiogenic TAMs vary between cancer types (SPP1+ in eight tumor types, VCAN+ in breast and melanoma, INHBA+ in stomach, FN1+ in kidney), so the function is conserved while the label is not; scoring continuous functional axes sidesteps this rather than colliding with it, and also avoids the multi-score assignment problem and the rare-category problem that had forced the exclusion of CD8.TPEX in Bloc6. Two non-standard symbols in the source table, IL23 and FASL, were mapped to IL23A and FASLG; coverage was complete (16/16, 34/34, 25/25, 4/4) in both histologies.

Per-patient median scores were compared between MPR and non-MPR using the Wilcoxon rank-sum test with Benjamini-Hochberg correction across the four axes, with the same ten-cell minimum used throughout the project.

As a post hoc follow-up in LUSC, where the M1 axis differed between response groups, the M1 signature was decomposed into four functional registers (inflammatory cytokines, IFN-induced chemokines, co-stimulation and antigen presentation, tryptophan catabolism). Each register was scored with UCell, and per-patient medians were compared between MPR and non-MPR with a Wilcoxon rank-sum test (Script 07), using the same ten-cell minimum. No multiple-testing correction was applied across the four registers, and the decomposition was not planned in advance.

#### Subtype markers (Wang et al.)

The four macrophage subtype markers reported in a treatment-naive NSCLC cohort (FABP4, FCN1, SPP1, SELENOP; Wang et al., 2022) were tested by assigning each cell to its highest-expressing marker. This is not a reproduction of the source method, which used unsupervised clustering rather than a threshold on four genes, and whose authors report that marker expression was not exclusive to a single cluster; the test asks only whether these markers structure the compartment here. Cells expressing none of the four were left unassigned. Patients with no cell assigned to a given marker were counted as zero percent rather than excluded, since exclusion would have removed the most extreme values and biased the comparison, the same pitfall encountered in Bloc5.
Patients with fewer than 10 assigned cells in total were excluded. Per-patient proportions were compared between MPR and non-MPR with a Wilcoxon rank-sum test per marker group, with Benjamini-Hochberg correction across the four groups.

#### Metabolic programs (Hallmark)

Four MSigDB Hallmark gene sets (glycolysis, hypoxia, oxidative phosphorylation and fatty acid metabolism) were scored with UCell to test a metabolic hypothesis. They were selected based on Yang et al. (2021, *Front Immunol*), who report that two immunosuppressive TAM subtypes in NSCLC differ by metabolism (SPP1+ macrophages glycolytic, CCL18+ macrophages fatty-acid-oxidative), and because the macrophage compartment of this cohort is dominated by SPP1 (Script 08). Selection was fixed before any computation; the remaining forty-six sets were not tested. No directional prediction was made, since Yang et al. studied a treatment-naive cohort and their work does not establish whether this axis tracks response after anti-PD-1 treatment. Per-patient medians were compared between MPR and non-MPR with a two-sided Wilcoxon rank-sum test, with Benjamini-Hochberg correction across the four programs. Spearman correlations between the four programs were computed at the patient level.

### Results

#### Functional scoring (Cheng et al.)

- **Only the M1 axis differs between response groups, and it survives correction in LUSC only.** In LUSC (86 patients: 28 MPR, 58 non-MPR), M1 is lower in responding tumors (median 0.060 versus 0.088; adjusted P = 0.019). In LUAD (43 patients: 9 MPR, 34 non-MPR), the difference runs in the same direction (0.092 versus 0.105) but does not reach significance after correction (adjusted P = 0.158). No other axis (M2, angiogenesis, phagocytosis) reaches significance in either histology.

*LUAD (43 patients: 9 MPR, 34 non-MPR)*

| Axis | Median MPR | Median non-MPR | p | p_adj |
|---|---|---|---|---|
| M1_UCell | 0.0918 | 0.105 | 0.081 | 0.158 |
| Angiogenesis_UCell | 0.0686 | 0.0947 | 0.086 | 0.158 |
| M2_UCell | 0.152 | 0.167 | 0.118 | 0.158 |
| Phagocytosis_UCell | 0.401 | 0.408 | 0.692 | 0.692 |

*LUSC (86 patients: 28 MPR, 58 non-MPR)*

| Axis | Median MPR | Median non-MPR | p | p_adj |
|---|---|---|---|---|
| M1_UCell | 0.0601 | 0.0883 | 0.0048 | 0.0192 |
| M2_UCell | 0.164 | 0.164 | 0.331 | 0.441 |
| Phagocytosis_UCell | 0.370 | 0.402 | 0.226 | 0.441 |
| Angiogenesis_UCell | 0.0902 | 0.0946 | 0.516 | 0.516 |

- The phagocytosis axis relies on a substantially smaller gene set (4 genes) than the other three axes (16 to 34 genes each), which limits how reliably it can be scored; this is a caveat on signature size, not a mechanistic claim.

- **M1 decomposition (LUSC, post hoc).** To understand what drives the lower M1 score in responding tumors, the signature was decomposed into its functional registers. The Cheng et al. M1 signature (16 genes) mixes inflammatory cytokines (TNF, IL1A, IL1B, IL6, IL23A), IFN-induced chemokines that recruit T cells (CXCL9, CXCL10, CXCL11, CCL5), co-stimulation and antigen presentation (CD86, CD40, CCR7, IRF1, IRF5), and tryptophan catabolism (IDO1, KYNU), which is immunosuppressive despite being IFN-inducible. A lower overall score could therefore reflect reduced inflammation or reduced IFN-gamma signalling, two different biological readings. Only the co-stimulation register differed between groups (median 0.098 in MPR versus 0.132 in non-MPR; unadjusted p = 0.012, four registers tested without correction); the other three registers scored zero for the median cell in both groups, meaning they were undetected rather than differentially expressed.

#### Subtype markers (Wang et al.)

- **No marker group differs significantly between response groups, in either histology.** Each cell was assigned to the marker with the highest expression among FABP4, FCN1, SPP1 and SELENOP (unassigned if all four were undetected). Patients with fewer than 10 assigned cells in total were excluded. Wilcoxon rank-sum test per marker group, with Benjamini-Hochberg correction across the four groups.

*LUAD (40 patients: 8 MPR, 32 non-MPR)*

| Marker | Median % MPR | Median % non-MPR | p | p_adj |
|---|---|---|---|---|
| SPP1 | 32.7 | 66.1 | 0.051 | 0.202 |
| SELENOP | 32.5 | 15.4 | 0.211 | 0.233 |
| FCN1 | 11.6 | 7.81 | 0.223 | 0.233 |
| FABP4 | 1.20 | 0.745 | 0.233 | 0.233 |

*LUSC (82 patients: 25 MPR, 57 non-MPR)*

| Marker | Median % MPR | Median % non-MPR | p | p_adj |
|---|---|---|---|---|
| FABP4 | 4.82 | 0 | 0.077 | 0.205 |
| SPP1 | 45.3 | 59.4 | 0.103 | 0.205 |
| FCN1 | 12.8 | 13.6 | 0.884 | 0.948 |
| SELENOP | 17.2 | 18.9 | 0.948 | 0.948 |

In both histologies, the SPP1-dominant fraction is numerically higher in non-MPR than in MPR tumors, but the difference does not survive correction in either.

#### Metabolic programs (Hallmark)

- **No program differs between response groups, in either histology.**

*LUAD (43 patients: 9 MPR, 34 non-MPR)*

| Program | Median MPR | Median non-MPR | p | p_adj |
|---|---|---|---|---|
| Fatty acid metabolism | 0.0781 | 0.0845 | 0.026 | 0.105 |
| Hypoxia | 0.131 | 0.139 | 0.133 | 0.267 |
| Glycolysis | 0.0722 | 0.0747 | 0.243 | 0.324 |
| Oxidative phosphorylation | 0.157 | 0.161 | 0.471 | 0.471 |

*LUSC (86 patients: 28 MPR, 58 non-MPR)*

| Program | Median MPR | Median non-MPR | p | p_adj |
|---|---|---|---|---|
| Glycolysis | 0.0759 | 0.0799 | 0.106 | 0.423 |
| Fatty acid metabolism | 0.0874 | 0.0872 | 0.291 | 0.583 |
| Hypoxia | 0.131 | 0.135 | 0.622 | 0.793 |
| Oxidative phosphorylation | 0.178 | 0.181 | 0.793 | 0.793 |

In both histologies, the small differences that appear follow no consistent pattern that a glycolytic versus oxidative dichotomy would predict. Fatty acid metabolism, glycolysis and oxidative phosphorylation are strongly and positively correlated at the patient level (Spearman rho 0.56 to 0.82), rather than behaving as opposing axes; only hypoxia is comparatively independent of the other three. These four programs therefore do not provide a usable functional readout in this cohort.

#### Cross-framework synthesis

- **Across the three frameworks, only the Cheng M1 axis differs between response groups after correction, and only in LUSC** (adjusted P = 0.019; LUAD, adjusted P = 0.158). The Wang subtype markers and the Hallmark metabolic programs show no significant difference in either histology.

- **Where a trend exists, it runs in the same direction in both histologies.** Cheng M1 and the SPP1-dominant fraction are numerically higher in non-MPR than in MPR tumors in both LUAD and LUSC, but only the LUSC M1 result survives correction. These observations are descriptive and are not interpreted further.