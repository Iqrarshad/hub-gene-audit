# Analysis history

> Scope note: this repository is the hub-gene centrality audit pipeline.
> Some early entries below reference exploratory stages that were not
> carried into this release. The stages present in `R/` are the audit
> pipeline, numbered 01 to 38 in dependency order.


Recorded so that withdrawn results are not cited by mistake, and so the
reasoning behind each correction is auditable.

The withdrawn scripts are not shipped. This file is the record of what they
did and why they were replaced. Stages 20, 32, 39, 40, 45, 46 and 49 have no
corresponding file in `R/`, and the gaps in stage numbering are deliberate.

## Literature analysis

| Version | Studies | Background | Outcome | Status |
|---|---|---|---|---|
| Stage 32 | 9 | 247 | p = 0.895 | Invalid. The null drew each study's list from a shared pool, which makes recurrence depend only on pool size. Sampling from any 25-gene pool gives the same expected recurrence whether the pool is high-degree or random. |
| Stage 40 | 9 | 247 | p = 0.112 | Superseded. Per-gene degree matching fixed the pool problem, but the background was too small. |
| Stage 45 | 30 | 247 | p = 0.0004 | Superseded. The random null returned 65.6 against an observed 44; a null larger than the observation cannot be interpreted. |
| Stage 46 | 30 | 1,274 | p_more = 1.0 | Superseded. Both nulls moved in opposite directions relative to stage 45, showing the earlier result was an artefact of background size. Gate L1 was set on requested rather than resolved genes. |
| **Stage 47** | **30** | **1,909** | **p_more = 0.949, p_fewer = 0.087** | **Report this.** All four gates pass; match error 3.7%. |

The direction effect was stable across every background size tested:
upregulated hub lists at the 96.2nd-96.8th degree percentile, downregulated
at 75.7th-78.1st, p ~ 1e-10. It is a property of the hub genes, not of the
null, and does not depend on the background.

## Corrections attempted, in order

| Stage | Approach | Bias AUC | Outcome |
|---|---|---|---|
| 19 | Binomial degree correction | 0.825 | Improves on 0.994 but far from 0.5 |
| 25 | Rewiring z-score, max | 0.779 | Cherry-picks each gene's extreme metric |
| 25 | Rewiring z-score, mean | 0.533 | Nominates only degree-1 and 2 genes; fails a degree floor |
| 33 | Composite, four rules | 0.729-0.829 | No better than the binomial alone |
| 34 | Reproducibility-weighted eigenvector | 0.593 | Reproducibility weighting beats plain correlation (0.715) but still recovers degree at r = 0.84 |
| 29 | Cohort co-expression network | 0.576 | Substrate change, the best of the scoring-side attempts |
| 42 | Citation-penalised score | 0.745 | Best variant keeps six of eight conventional hubs |
| 43 | CRISPR co-essentiality network | 0.483 | Bias eliminated, but selectivity 0.552, i.e. chance |
| 35 | Full redesigned pipeline | 0.486 | Zero enriched GO terms |
| 44 | Propagation gates | - | Kill criterion met: RWR reduces citation correlation rather than amplifying it |

## Validations that failed, and why

**GO enrichment (stage 35).** Circular. GO annotation exists for genes that
have been studied, so a method designed to escape annotation density selects
genes that by construction have no annotations. The test cannot distinguish
"found nothing" from "found uncharacterised genes".

**Survival (stages 36, 37).** Direction reversed across cohorts within IDH
strata after age and grade adjustment (HR 0.697, 1.725, 1.353). Size-matched
random gene sets reached p < 0.05 in 31-52% of draws, so the significance
was uninformative.

**DepMap selectivity (stage 38).** No gene set separated better than chance
(AUC 0.510-0.572, all p > 0.07), so it cannot serve as a calibration target.

## Errors corrected during development

- Stage 06 and 17 computed `stress` and `radiality` as duplicates of
  betweenness and closeness, inflating the consensus count by two free votes
- Stage 09 initially adjusted for purity only; the immune covariate was
  computed and dropped from the design matrix
- Stage 28 established that the split-half stability ceiling is 0.51, not
  the 0.33 estimated by hand; earlier instability conclusions were withdrawn
- Stage 38's DOS viability test ran on 148 genes rather than the
  transcriptome, so its R^2 was range-restricted
- Stage 44's first mapping used the STRING aliases file, which carries many
  alias types per protein and reduced the usable set to 65 genes;
  `protein.info` gives one preferred name per protein

## File and naming changes

Renames made when the archive was prepared for deposit. Stage numbers,
function names and analysis content are unchanged.

| Stage | Old file | New file |
|---|---|---|
| 18 | `18_enrichment_null_and_hla.R` | `18_connectivity_null_and_hla.R` |
| 38 | `38_dos_and_depmap.R` | `38_depmap_validation.R` |
| 44 | `44_apnp_gate0.R` | `44_propagation_gates.R` |
| 47 | `47_literature_final_v3.R` | `47_literature_scale_analysis.R` |

Stage 44 outputs were renamed with the file: `apnp_gate0.csv` and
`apnp_gate0.rds` are now `propagation_gates.csv` and
`propagation_gates.rds`. The APNP acronym was dropped throughout in favour
of the written form. Results from an earlier run are not readable by stage
number alone and should be regenerated.

Stage 48 and the `figures` group were added after the first packaged
version, which stopped at stage 47.

## Other corrections at packaging

- Stage 41 treated a malformed stage 47 degree table the same as a missing
  one, dropping genes from the decomposition without a message. Absence is
  now logged; a file present without a degree column is an error.
- Superseded stage 39 was described in the runner and README as a channel
  comparison. It tests recurrence on the physical-evidence-only network.

## Propagation audit

Stage 44 reported rho(x*, E) = 0.855 and this was read as retention of the
expression signal. The reading is close to circular: x* = A^-1 E contains E
by construction, so the correlation largely confirms that diffusion did not
destroy its own restart vector. Stages 50 and 51 replace it with held-out
evaluation.

| Stage | Question | Result |
|---|---|---|
| 49 | Does a citation-orthogonal operator help? | Withdrawn. Criteria had no benefit term, so the sweep selected alpha near zero, which is no diffusion at all |
| 50 | Why does propagation recover degree? | Held-out recovery 0.065 against true values, 0.71 against degree. Degree predicts top-decile effect size at AUC 0.417, below chance; neighbourhood mean at 0.649 |
| 50 | Is the signal graph-structured at all? | Autocorrelation 0.217 against a degree-stratified null of 0.003, z = 19.7 |
| 50 | Is 0.065 low in absolute terms? | Calibration curve over signals of graded smoothness puts expected recovery at 0.259. Degree-adjusted recovery is 0.198, or 76% of that |
| 51 | Does constraining the operator beat adjusting the output? | No. Best constrained 0.158 against 0.198 post-hoc, at every alpha. Rank association with degree stays at 0.41 with the constraint satisfied exactly |

A single linear constraint removes the mean projection onto the degree
direction and leaves the rank association intact, which is the quantity
that matters. The optimal penalty is the projection limit, so there is no
tuning parameter and no method. Post-hoc degree adjustment is reported.

Three thresholds were set and all three were wrong on first use: K5 in
stage 49 had no benefit term, D6 in stage 50 used an absolute cut of 0.5
when the control family could not exceed 0.36, and the first version of D6
built a degenerate control from 200 sparse seeds. Absolute thresholds on
quantities whose achievable range has not been computed are unreliable;
the calibration curve replaced the last of them.

## Path configuration

`DATA_DIR` and `RESULTS_DIR` were hardcoded to absolute Windows paths on one
machine, which made the repository unusable elsewhere and exposed a local
directory layout. They are now resolved in order from a command line flag,
an environment variable, an untracked `config.local.R`, then `data/` and
`results/` beside the project root. The pipeline stops if the two resolve to
the same directory and warns if `DATA_DIR` does not exist. Resolved paths
are printed at the start of every run.

`run_all.R` gained `--clean`, which empties the derived output directories
and reruns every stage. `DATA_DIR` is not touched, so downloads are not
repeated.

## Inputs

Source datasets, none redistributed here.

| Dataset | Accession | Role |
|---|---|---|
| Glioma expression, RNA-seq | GSE147352 | discovery |
| Gravendeel microarray | GSE16011 | discovery |
| REMBRANDT | GSE108474 | discovery |
| Glioma microarray | GSE15824 | discovery |
| Glioma microarray | GSE21354 | discovery |
| TCGA-LGG | cBioPortal / GDC | validation |
| TCGA-GBM | cBioPortal / GDC | validation |
| CGGA mRNAseq 325 | CGGA.org.cn | validation |
| CGGA mRNAseq 693 | CGGA.org.cn | validation |
| Single-cell glioma | GSE84465 | ground-truth reference |
| STRING human network | STRING v12 | interactome |
| Gene-publication links | NCBI gene2pubmed | citation counts |
| CRISPR essentiality | DepMap | co-essentiality substrate |

## Path configuration

DATA_DIR and RESULTS_DIR are resolved from, in order, a command-line flag,
an environment variable, an untracked config.local.R, then data/ and
results/ beside the project root. The pipeline stops if the two resolve to
the same directory and warns if DATA_DIR is absent. Resolved paths are
printed at the start of every run. run_all.R accepts --clean to empty the
derived output directories and rerun; downloaded source data is left in
place.

## Added stages 42 to 44 (interval, robustness, corrections)

These stages are additive. They read the cached outputs of 01 to 41 and
write new tables only; no earlier stage or output is modified, so a rerun
reproduces the existing results unchanged and adds the new files.

- **42 bootstrap_intervals** writes `bootstrap_intervals.csv`: BCa intervals
  with 5000 resamples for the headline AUCs that were previously reported as
  point estimates only (publication and residual prediction of hub
  membership, the 131 published-hub AUCs, the enrichment analogue, and the
  DepMap glioma-selective dependency of the hub set). The candidate-network
  degree-to-hub AUC of 0.994 keeps its stage 12 DeLong interval and is not
  re-bootstrapped. AUC is the tie-corrected Mann-Whitney statistic; the BCa
  jackknife is computed in closed form from placement values.
- **43 multiseed_robustness** writes `multiseed_robustness.csv` and
  `multiseed_robustness_raw.csv`: each interval recomputed across seeds
  (default 20) to show the endpoints are stable to the resampling seed. The
  point estimates are deterministic and do not move with the seed. Seed
  stability of the STRING spike-in null (stage 07) is separate and is not
  triggered here, because re-drawing its backgrounds calls the STRING API.
- **44 corrections_table** writes `corrections_all.csv`: the nine corrections
  in two families plus the two network-propagation variants (random walk
  with restart, degree-constrained diffusion) and the cytoHubba baseline, in
  one table. Reproducibility and attention-dependence apply to the nine;
  propagation is judged on recovery of held-out signal versus degree, so
  those cells are native to propagation and the shared columns are left NA.
  Every numeric cell is read from a stage output; a value that cannot be
  found is left NA.

Helper `00h_bootstrap.R` holds the shared base-R functions (`fast_auc`,
`bca_auc_ci`, tolerant readers). Base R only, no new package dependency.
Configurable constants: `N_BOOT_CI` (stage 42), `MS_SEEDS` and `MS_R`
(stage 43).

## Reproducibility revision (Methods-to-code reconciliation)

Applied to make the code match the Methods as written, ahead of the final run:
- 00_config.R: deg_min_datasets kept at 3. A strict five-way intersection (=5) was
  tested and yields 0 composition-robust genes, so it is not viable; the pre-declared
  >=3-of-5 direction-consistent rule (148 genes) stands. The manuscript should describe it
  as "differentially expressed in at least three of the five discovery cohorts with
  consistent direction" rather than an "intersection".
- 17_coexpression_network.R: co-expression edges now taken on positive Spearman
  correlation (was |rho|); log message updated.
- 07_hub_null.R: spike-in null now uses the 6-metric, >=4 consensus rule that the
  main pipeline uses (removed the duplicate stress/radiality metrics and the >=6
  threshold).
- 02_deg_microarray.R: REMBRANDT grade contrast (GSE108474_HGG_vs_LGG) now adjusts
  for neuronal + glial + immune, matching Eq. 5 and stage 41.
- 24_literature_bias_gene2pubmed.R: degree->publication regression now retains
  zero-degree and zero-publication genes (left join, zero-fill), matching the
  stated log(1+.) rationale.
- 36_wgcna_reproducibility.R: writes per-cohort soft-thresholding power (beta) to
  tables/wgcna_soft_threshold_power.csv.

## STRING density-cache fix
- 14_density_threshold.R: string_net() no longer caches a NULL/empty result
  (a rate-limited STRING call was poisoning cache/string_density/t<score>.rds).
- 15, 17, 18, 25: if the density cache is NULL/empty, fall back to the stage-06
  edge cache rds/string_<score>.rds. Prevents "No cached STRING edges" failures.

## STRING rate-limit robustness (stage 14)
- 14_density_threshold.R: string_net() now reuses the stage-06 edge cache
  (rds/string_<score>.rds) for scores already fetched, and retries with
  exponential backoff (5-60s) otherwise. Fixes the confidence-threshold sweep
  returning 0 edges at every score when STRING rate-limits after stages 06-13.
- 08_centrality_diagnostics.R: log corrected "eight" -> "six" centrality metrics.

## Stale-cache self-healing (string_density)
- --clean clears rds/tables/figures but not cache/, so a NULL string_density
  cache written by a pre-fix run persisted and was read first.
- 14_density_threshold.R: string_net() now ignores a NULL/empty cached file
  and refetches (or reuses the stage-06 cache).
- 41_reproducibility_gated_selection.R: falls back to rds/string_<score>.rds
  when the density cache is NULL/empty (fixes the per-cohort conventional row
  reading 0 genes / NA).

## Repo finalization (pre-Zenodo hygiene)
- README: removed TCGA-GBM (not used); corrected 54.5%/28.4% -> 47.9%/30.7%
  (and "text-mining and co-occurrence channels"); "Eight" -> "Nine" corrections
  with the correct list; figure stages 48/52/56 -> 29/32/38/41.
- 44_corrections_table.R: binomial specific-connectivity dependence AUC now
  read from hubbias_predictability.csv (was NA / "(text)"); matches Table 2.
- published_hub_lists.csv: added cancer_group column (21 cancer strings collapse
  to the 15 organ-site groups the manuscript reports).
- make_supplementary_tables.R: new consolidation script that maps outputs to the
  manuscript S1-S10 supplementary-table names. Run after run_all.R.
