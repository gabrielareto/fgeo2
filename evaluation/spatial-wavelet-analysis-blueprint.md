# Spatial wavelet analysis: deferred development blueprint

Date: 2026-09-01  
Status: design record and future roadmap; **not the active MVP scope**  
Legacy target: `03_spatial/R/new/spatial-block.analysis.r`  
Evaluation basis: repository commit `040b73e0cd9a0d55d64c1639b22159d5c6cbd20e`

## 1. Purpose and scope decision

Spatial wavelets could eventually provide a coherent family of scale-specific
analyses for ForestGEO and the broader ATFS community. The family could cover
one or many species, repeated censuses, environmental fields, individual marks,
and two- or three-dimensional domains.

That is a useful architectural horizon, but it is too broad for the first
deliverable. This document preserves the horizon so that a small MVP is named,
structured, and documented in a way that can be extended without redesigning
its vocabulary or mixing incompatible statistics.

The roadmap is therefore governed by two rules:

1. implement and validate the smallest scientifically defensible analysis first;
2. treat every later capability as a separately justified increment, not as an
   implied commitment of the MVP.

## 2. Intended users and design consequences

The intended users are ForestGEO researchers, technicians and data managers,
Smithsonian staff, and researchers in the wider ATFS community. The code should
help them answer recurring ecological questions while remaining usable outside
ecology.

This implies:

- general statistical function names, with ecological use cases in tutorials;
- transparent calculations and ordinary R objects that can be inspected;
- small functions with one responsibility each;
- minimal, focused, seriously maintained dependencies;
- explicit scientific assumptions and unsupported cases;
- deterministic examples and validation scripts;
- null-pattern generation separated from calculation of the statistic;
- a common vocabulary across input types and ecological examples;
- no species-, census-, or site-specific wrapper when `wavelet_by_group()` can
  express the same workflow.

## 3. Active MVP boundary

The first implementation should evaluate only a **global, isotropic, two-
dimensional wavelet variance for one unmarked point pattern**.

The MVP candidate is bounded as follows:

| Dimension | MVP decision |
|---|---|
| Input | One `unmarked_point` pattern |
| Domain | Complete rectangular 2-D window |
| Statistic | Global normalized wavelet variance |
| Transform | Isotropic Morlet filtering in the Fourier domain |
| Boundary | Periodic/toroidal calculation, stated explicitly |
| Baseline | Homogeneous CSR expectation of one |
| Inference | Corrected analytical CSR interval only if independently validated; otherwise external simulations |
| Null generation | Outside the statistic function |
| Output | Ordinary data frame, plus a small metadata list if needed |
| Tutorial | One ecological example, one null-envelope example, and direct comparison with familiar spatial summaries |
| Explicitly excluded | Multiple groups, bivariate statistics, marked points, fields, local maps, anisotropy, irregular windows, gaps, and 3-D |

The MVP should not expose a broad public API before its normalization, scale
calibration, boundary behavior, and inference have passed the validation gates
in section 11.

## 4. Shared vocabulary

Use the following terms consistently in code, documentation, tests, and issue
discussion:

- **wavelet spectrum**: scale-indexed spectral quantity before a more specific
  summary is named;
- **wavelet variance**: univariate scale-specific variation;
- **wavelet covariance**: signed cross-variation retaining physical units;
- **wavelet correlation**: normalized signed association in `[-1, 1]`;
- **wavelet coherence**: normalized magnitude-squared association in `[0, 1]`;
- **pattern**: a spatial point pattern;
- **unmarked_point**: locations without a quantitative individual mark;
- **marked_point**: locations with a quantitative or categorical individual mark;
- **field**: a value observed throughout a spatial grid or continuous domain;
- **mark**: an attribute attached to an individual point;
- **envelope**: a summary calculated from already-generated null statistics.

Do not use `coherence` for the signed, unbounded transformation returned by the
legacy `wavelet.bivariate()` function. Do not use `points` as shorthand for
`unmarked_point` when the distinction from `marked_point` matters.

## 5. Systematic analysis map

The table is an architectural map, not an implementation checklist. Each row
requires its own ecological hypothesis, normalization, null hypothesis,
boundary analysis, and validation before it becomes supported.

| Ecological question | Input A | Input B | Primary summary | Typical grouping or contrast | Candidate null hypothesis |
|---|---|---|---|---|---|
| At what scales is one species aggregated or overdispersed? | `unmarked_point` | — | variance | one species | homogeneous or fitted inhomogeneous CSR |
| Do aggregation scales differ among species? | `unmarked_point` | — | variance | species | group-specific null patterns |
| Do aggregation scales differ among censuses? | `unmarked_point` | — | variance | census within species | census-specific null patterns |
| Do life stages or size classes occupy different scales? | `unmarked_point` | — | variance | size class within species | random labelling or class-specific process |
| Does a population retain spatial structure through time? | `unmarked_point` | `unmarked_point` | covariance or correlation | census pair within species | temporal randomization chosen for the process |
| Are two species associated at particular scales? | `unmarked_point` | `unmarked_point` | covariance or correlation | species pair | independence, random labelling, or fitted-process null |
| Is community-scale structure shared among many species? | `unmarked_point` | `unmarked_point` | pairwise or aggregate cross-summary | species pairs or guilds | explicit community null model |
| Does abundance track an environmental gradient by scale? | `unmarked_point` | `field` | cross-spectrum, covariance, or correlation | species × environmental variable | point-field independence preserving relevant marginals |
| Do species respond similarly to an environmental field? | `unmarked_point` | `field` | cross-summary | species | group-specific point-field nulls |
| At what scales is individual size or biomass structured? | `marked_point` | — | marked variance | species, census, or size class | random labelling or a fitted marked process |
| Do marks persist through censuses? | `marked_point` | `marked_point` | marked covariance or correlation | census pair | mark reassignment or temporal null |
| Are individual marks associated with environment? | `marked_point` | `field` | marked point-field cross-summary | species or census | random labelling and/or point-field independence |
| Do two spatial fields share scale-specific structure? | `field` | `field` | covariance, correlation, or coherence | field pair or time pair | surrogate fields preserving stated autocorrelation |
| Does a field retain structure through time? | `field` | `field` | covariance, correlation, or coherence | time pair | temporal or surrogate-field null |
| Is spatial structure directional? | any supported type | optional second input | directional spectrum or cross-spectrum | direction | rotation-appropriate null |
| Where, rather than only at what scale, does structure occur? | any supported type | optional second input | local wavelet map | location and scale | location-aware simulation envelope |
| Does vertical structure add information beyond mapped x-y position? | 3-D pattern or field | optional second input | dimension-general summary | layer, height, or time | 3-D process-specific null |

The null hypotheses in the last column are deliberately not API defaults. For
many rows, several nulls answer different ecological questions.

## 6. Candidate function family

Names below reserve a coherent vocabulary. Only functions needed by a validated
increment should be implemented. Avoid placeholder exports.

### 6.1 Input preparation

| Candidate function | Responsibility | Earliest increment |
|---|---|---|
| `rasterize_unmarked_point()` | Convert an unmarked point pattern to cell counts while retaining window and cell metadata | MVP, internal first |
| `rasterize_marked_point()` | Convert a marked point pattern to explicit mark summaries per cell | marked-point increment |
| `prepare_wavelet_field()` | Validate and standardize a gridded field without silently imputing gaps | field increment |

These functions must not hide the distinction between counts, mark sums, mark
means, or other cell summaries. A marked compound point process does not share
the unmarked CSR normalization: for compound Poisson cell sums,
`Var(X) / E(X) = E(Z^2) / E(Z)`, not one.

### 6.2 Numerical primitives

| Candidate function | Responsibility | Earliest increment |
|---|---|---|
| `wavelet_scales()` | Construct and report the evaluated scale grid and Fourier-period conversion | MVP |
| `wavelet_filter()` | Construct the selected wavelet filter for a grid and scale | MVP, possibly internal |
| `wavelet_spectrum()` | Calculate a univariate scale-indexed spectrum | MVP |
| `wavelet_cross_spectrum()` | Calculate a scale-indexed cross-spectrum | first bivariate increment |

The numerical primitives should expose enough intermediate values for an
independent equation-to-code audit. They should not choose an ecological null.

### 6.3 Statistical summaries

| Candidate function | Responsibility | Input combinations |
|---|---|---|
| `wavelet_variance()` | Calculate a univariate scale-specific variance | one supported input |
| `wavelet_covariance()` | Calculate signed cross-variation in physical units | two supported inputs |
| `wavelet_correlation()` | Calculate normalized signed association | two supported inputs |
| `wavelet_coherence()` | Calculate magnitude-squared association | two supported inputs |
| `wavelet_map()` | Calculate a local, location-by-scale result | later local-analysis increment |
| `wavelet_anisotropy()` | Calculate a direction-specific result | later anisotropy increment |

The same statistical names should eventually work across supported input
combinations. Input-type-specific normalization and validation should be
explicit; it must not be inferred by a fragile heuristic.

### 6.4 Collections, inference, and display

| Candidate function | Responsibility | Constraint |
|---|---|---|
| `wavelet_by_group()` | Apply one supported statistic independently by species, census, size class, guild, site, or another grouping variable | tutorial supplies ecological groupings; no `wavelet_by_species()` wrappers |
| `wavelet_interval_csr()` | Calculate only a validated analytical homogeneous-CSR interval | not a null generator |
| `wavelet_envelope()` | Summarize observed and already-calculated null statistics | accepts null results; never simulates patterns |
| `plot_wavelet()` | Display common result columns without changing the statistic | thin convenience layer only after result schema stabilizes |

Null patterns or fields should be generated by maintained, domain-appropriate
software or by user code, then passed through the same statistic pipeline as
the observed data. fgeo2 should not hide null generation inside
`wavelet_variance()`, `wavelet_correlation()`, or `wavelet_envelope()`.

## 7. Common result contract

Where possible, every global statistic should return an ordinary data frame
with stable columns such as:

- `scale` and its coordinate unit;
- `statistic` and the precise statistic name;
- `value`;
- `input_type_a` and, when relevant, `input_type_b`;
- group identifiers supplied by `wavelet_by_group()`;
- explicit boundary and normalization identifiers;
- interval or envelope columns only when calculated by a separate inference
  function.

Small metadata that cannot be represented safely per row may be returned in a
plain list with the data frame. A custom class is justified only if it clearly
improves validation or interoperability. Raw spectra, filters, and other audit
objects should be optionally recoverable without making every ordinary result
large.

## 8. Dependency and software-resource audit

No dependency should be adopted merely because it contains the word
"wavelet". Each candidate receives one of four roles:

1. required dependency;
2. optional infrastructure;
3. validation comparator;
4. rejected, with a recorded reason.

The audit should record, at minimum:

- package and version inspected;
- current maintainers and authoritative source;
- release recency and continuity;
- issue and pull-request activity;
- automated checks and supported R versions;
- license compatibility;
- exact functions relevant to this work;
- scientific definition of those functions;
- whether the package is required at run time or only during validation;
- reproducibility and numerical limitations.

Starting candidates include base R's FFT, the maintained `spatstat` family for
point-pattern representation and simulation, and `waveslim` as a comparator for
some raster-field questions. Raster infrastructure should be selected only
after comparing maintained candidates against the actual field requirements.
Archived, abandoned, or opaque code may inform provenance but should not become
a production dependency.

## 9. Provenance and credit, including Matteo Detto's resources

Before reimplementing beyond the MVP, search systematically for improved or
authoritative code through:

- the publisher's supplementary material for the method papers;
- Smithsonian/ForestGEO publication and data repositories;
- Matteo Detto's institutional page, public website, and public code profiles;
- Tania Brenes's attributable resources;
- the public `forestgeo/ctfs` repository and commit history;
- later papers that reused or modified the calculations;
- archived package releases and source bundles.

For every discovered resource, record its URL, retrieval date, authorship,
license, relationship to the legacy code, and whether equations and outputs can
be reproduced. If no corrected research implementation is public, prepare a
short author query describing the exact uncertainties. Do not send it without
the repository owner's authorization.

Credit belongs in code comments where an equation is implemented, in function
documentation, in the tutorial, and in a bibliographic file. Reimplementation
must not erase the provenance of the Detto, Muller-Landau, Brenes, and CTFS work.

## 10. Literature review plan

The literature review has two linked tracks.

### 10.1 Statistical method

Record the definitions and assumptions of:

- global and local wavelet spectra;
- isotropic and directional transforms;
- point-process normalization;
- cross-spectrum, covariance, correlation, and coherence;
- effective degrees of freedom and analytical intervals;
- finite grids, padding, masks, gaps, and periodic boundaries;
- surrogate-field methods and simulation envelopes;
- marked point processes and compound-Poisson baselines;
- extensions from 2-D to 3-D.

### 10.2 Ecological questions

Record examples involving:

- aggregation and scale separation;
- interspecific association;
- temporal formation, persistence, and dissipation of clumping;
- topographic and environmental association;
- canopy or remotely sensed fields;
- individual size, biomass, traits, or demographic status;
- habitat association versus dispersal or density dependence;
- multi-species and community comparisons.

Each source should be classified as theoretical, empirical, software,
interpretive, or unresolved. The review should state what a result can and
cannot establish ecologically; a scale-specific association is not by itself a
causal mechanism.

Initial primary anchors are:

- [Detto & Muller-Landau (2013), *Fitting Ecological Process Models to Spatial Patterns Using Scalewise Variances and Moment Equations*](https://www.journals.uchicago.edu/doi/full/10.1086/669678);
- [Detto et al. (2013), *Tropical Forest Canopy Disturbances and Recovery Are Driven by Competition for Light and Proximity to Streams*](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0076296);
- the public [`forestgeo/ctfs` implementation history](https://github.com/forestgeo/ctfs/commits/master/R/block.analysis.R).

These are starting points, not a complete review.

## 11. Validation and audit plan

Every implementation increment requires a deterministic validation harness
independent of its tutorial. Fixed seeds, dependency versions, code commit, and
machine-readable results must be recorded.

### 11.1 Equation-to-code checks

- map each implemented equation to the exact code expression and unit test;
- calculate simple filters and spectra independently when feasible;
- test Parseval/scaling identities relevant to the chosen normalization;
- expose the squared-filter effective-degrees-of-freedom calculation;
- test invariants such as circular-translation invariance under a periodic FFT;
- distinguish mathematical scale, Fourier period, grid cells, and coordinate units.

### 11.2 Simulation matrix

The MVP validation should cross representative values of:

- fixed-abundance CSR and Poisson CSR;
- low, medium, and high abundance;
- square and rectangular grids;
- odd and even grid dimensions;
- at least two cell sizes per window;
- complete scale ranges and conservative truncated ranges;
- clustered, inhibited, two-scale, and inhomogeneous processes;
- patterns translated toward or across boundaries;
- known regression cases from the legacy function.

Later increments add independent, positively associated, negatively associated,
identical, lagged, marked, point-field, field-field, masked, anisotropic, and
3-D cases as appropriate.

### 11.3 Behavioral expectations

- normalized unmarked-point variance is centered near one under the stated CSR
  convention;
- empirical interval coverage matches its stated level within Monte Carlo error;
- known cluster and inhibition processes depart in the expected direction and
  scale range;
- bivariate independence is centered at the statistic's correct null value;
- correlation remains in `[-1, 1]` and coherence in `[0, 1]`;
- identical nondegenerate inputs reach their mathematically expected limit
  without undocumented `Inf` or `NaN`;
- observed and null data use the identical rasterization, transform, scale, and
  boundary pipeline;
- unsupported masks, units, or dimensional combinations fail explicitly.

### 11.4 Independent audit artifacts

When implementation begins, maintain:

- a rendered evaluation QMD describing claims, equations, and evidence;
- a standalone validation R script runnable without the tutorial;
- a source/provenance ledger;
- a package audit with accepted and rejected candidates;
- machine-readable validation results;
- an equation–code–test crosswalk;
- a tutorial QMD with deterministic ecological cases;
- `sessionInfo()` or an equivalent dependency record;
- small commits and a pull request whose claims link to the artifacts.

Another human or AI should be able to audit a claim without relying on chat
history, hidden state, or an unrecorded manual step.

## 12. Incremental roadmap and gates

| Increment | Deliverable | Entry gate | Exit gate |
|---|---|---|---|
| R0 — provenance | Source ledger, literature map, package audit | blueprint accepted | method and software uncertainties are explicit |
| M0 — MVP | 2-D global variance for one `unmarked_point` pattern | R0 complete | equations, simulations, API, tutorial, and independent rerun pass |
| M1 — collections | `wavelet_by_group()` over the validated M0 statistic | stable M0 result contract | multiple species/census/class examples need no special wrappers |
| M2 — unmarked bivariate | covariance and/or correlation for two unmarked patterns | ecological null and target statistic selected | range, null behavior, and pair simulations pass |
| M3 — fields | field variance and selected point-field/field-field cross-statistics | maintained raster infrastructure selected | units, masks, surrogates, and boundary behavior pass |
| M4 — marked points | marked variance and selected cross-statistics | mark estimand and baseline derived | compound-process and random-labelling tests pass |
| M5 — local/directional | local maps and/or anisotropy | demonstrated recurrent use case | localization, multiplicity, direction, and edge tests pass |
| M6 — 3-D | dimension-general selected statistics | real 3-D ForestGEO/ATFS use case and data representation identified | scale, geometry, memory, visualization, and null validation pass |

An increment may be stopped, deferred, or rejected at its entry gate. Completing
M0 does not obligate M1–M6.

## 13. Tutorial architecture

The MVP tutorial should stay short and answer one ecological question well. As
increments are accepted, the tutorial can grow by reusable cases rather than by
function-by-function exposition:

1. one species: aggregation scales;
2. many species with `wavelet_by_group()`;
3. one species through censuses;
4. two species;
5. one species versus an environmental field;
6. marked individuals;
7. two fields;
8. constructing null patterns or fields externally and calculating an envelope;
9. interpreting agreement and disagreement with K, pair-correlation, quadrat,
   or other familiar summaries;
10. boundaries, scale limits, unsupported inputs, and causal caveats.

Every case should state the estimand, input type, null hypothesis, units,
boundary assumption, expected output, and ecological interpretation.

## 14. Stop conditions

Do not expand the implementation when any of the following holds:

- the ecological question is not more clearly answered than with established
  maintained methods;
- the target statistic or normalization is ambiguous;
- relevant author code or a maintained implementation already solves the need
  more reliably;
- edge, mask, or scale behavior cannot be bounded transparently;
- the only available dependency is unmaintained or scientifically opaque;
- validation requires an unverifiable reference result;
- the public API would imply unsupported combinations;
- tutorial complexity exceeds the practical benefit for intended users.

In that case, document how the current primitive could be extended, the
principles the extension must preserve, and the evidence still required. Leave
the implementation to future contributors rather than shipping an uncertain
generalization.

## 15. Immediate next step

Return to the MVP-first evaluation. The next work product should be the small R
validation harness described in the Stage 1 evaluation: one experimental
unmarked-point variance kernel, transparent filters and degrees of freedom,
fixed-seed simulations, and compact diagnostic outputs. No bivariate, marked,
field, grouped, local, anisotropic, or 3-D production function should be added
until the M0 exit gate is met.
