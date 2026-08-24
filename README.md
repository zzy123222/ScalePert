# ScalePert

[![tests](https://github.com/zzy123222/ScalePert/actions/workflows/tests.yml/badge.svg)](https://github.com/zzy123222/ScalePert/actions/workflows/tests.yml)

**ScalePert: interpretable local-manifold virtual perturbation for target prioritization from observational single-cell transcriptomes**

ScalePert is an interpretable, multiscale framework that uses naturally occurring transcriptional heterogeneity in observational single-cell (or single-nucleus) RNA-seq data to infer target-reduction-associated cellular-state transitions and prioritize candidate targets — without requiring perturbation-trained models, inferred regulatory networks, or additional molecular modalities.

The package implements the two modules described in the manuscript:

- **ScalePert-Cell** — embeds cells in a shared PCA space, builds population-restricted kNN graphs, shifts each cell toward transcriptionally similar neighbors with lower observed expression of a nominated target (averaging displacement vectors across a target set *before* functional reconstruction for multi-target hypotheses), reconstructs functional program scores at displaced coordinates from nearby observed cells, and quantifies the response with signed Wasserstein-1 distribution shifts.
- **ScalePert-Tissue** — maps population-level ScalePert-Cell outputs onto a ligand–receptor-derived communication graph and applies matrix-exponential (`exp(−βL)`) propagation to generate multicellular tissue-level priorities.

## Installation

```bash
pip install git+https://github.com/zzy123222/ScalePert.git
```

Or from a local clone:

```bash
git clone https://github.com/zzy123222/ScalePert.git
cd ScalePert
pip install .
```

Requirements: Python ≥ 3.9 with `numpy`, `pandas`, `scipy`, `scikit-learn`, `anndata`, `scanpy`, `statsmodels` (installed automatically).

## Quick start

### 1. End-to-end pipeline

```python
import scanpy as sc
from scalepert import ScalePertPipeline

adata = sc.read_h5ad("your_atlas.h5ad")          # requires .obs["cell_type"]
pipeline = ScalePertPipeline(cell_type_key="cell_type", k=30, scale=1.0, beta=0.3)
pipeline.fit(adata, targets=["TLR2", "ITGAM", "CCL5", "HIF1A", "VCAM1",
                             "JUN", "AIM2", "IFI16", "PYCARD"])

print(pipeline.cell_ranking())      # cell layer: most negative mean shift = highest reversal priority
print(pipeline.tissue_ranking())    # tissue layer: largest non-negative suppression score = top priority
pipeline.export("results/scalepert")
```

If your `AnnData` is already normalized, log-transformed and has `.obsm["X_pca"]` plus program scores in `.obs`, call `pipeline.fit(adata, score_existing=True)` to skip preprocessing.

### 2. Module-level usage

```python
from scalepert import (
    prepare_adata, score_programs, scalepert_cell,
    build_communication_graph, propagate_tissue,
)

adata = prepare_adata(raw_adata)                     # normalize → HVG → scale → PCA
scored = score_programs(adata)                       # built-in PROGRAMS by default,
                                                     # or pass {"my_program": [...gene list...]}
cell_res = scalepert_cell(scored, targets=["IFI16", "HIF1A"],
                          cell_type_key="cell_type", k=30, scale=1.0)
print(cell_res.ranking())
print(cell_res.table)                                # per-population signed W1 shifts

comm = build_communication_graph(scored, "cell_type")   # ligand–receptor graph
tissue_res = propagate_tissue(cell_res, comm, beta=0.3) # exp(−βL) propagation
print(tissue_res.ranking())
```

### 3. Example dataset

```python
from scalepert import load_example_data, ScalePertPipeline

adata = load_example_data()   # bundled synthetic kidney-style atlas
pipe = ScalePertPipeline().fit(adata)
print(pipe.cell_ranking())
```

## Score interpretation (read before ranking targets)

The two layers use **different score definitions and directions** and must be ranked independently; their numerical values are **not comparable across layers**.

| Layer | Quantity | Direction | Priority rule |
|---|---|---|---|
| ScalePert-Cell | signed mean Wasserstein-1 shift over the four programs | negative = disease-program suppression | more negative = higher priority |
| ScalePert-Tissue | propagated program-suppression score (total magnitude of downward-shifted programs per population, summed over populations) | non-negative by construction | larger = higher priority |

In disease-oriented applications where the predefined programs represent pathological processes, more-negative ScalePert-Cell shifts indicate stronger predicted suppression of disease-associated programs and therefore greater reversal priority.

## Method summary

Given a processed gene-by-cell matrix \(X\), population labels \(c(i)\), nominated targets and functional readouts:

1. **Shared representation** — normalize, log-transform, select highly variable genes, scale and compute PCA (50 components by default).
2. **Population-restricted neighborhoods** — cosine-distance kNN graph (k = 30 by default), constructed separately within each annotated population.
3. **Target-conditioned displacement** — for target gene *g*, the lower-expression neighbor set of cell *i* is
   \(\mathcal{N}_i^- = \{j \in \mathrm{KNN}(i) : e_j < e_i\}\),
   and the local displacement is
   \(\Delta_i = s \cdot \frac{1}{|\mathcal{N}_i^-|} \sum_{j \in \mathcal{N}_i^-} (x_j - x_i)\),
   with scaling factor \(s = 1.0\) (\(s = 1/q\) per target within a q-target combination). Cells without eligible lower-expression neighbors remain at their original coordinates.
4. **Additive combinations** — for a set G of q targets, individual displacement vectors are averaged **before** functional reconstruction.
5. **Functional reconstruction** — program scores at virtual coordinates are reconstructed as the unweighted mean of nearest observed cells in PCA space; cells displaced above the within-analysis 25th percentile are retained for distribution comparison.
6. **Signed Wasserstein-1 shift** — for program *m*:
   \(\mathrm{Shift}_m = \mathrm{sign}\!\left(\bar{Q}_m - \bar{P}_m\right) \cdot W_1(P_m, Q_m)\).
   Negative values denote suppression of disease-associated programs. ScalePert-Cell summarizes the four program-specific shifts by their arithmetic mean.
7. **Tissue aggregation** — annotated populations form nodes; candidate ligand–receptor interactions form directed edges weighted by Σ mean-ligand(sender) · mean-receptor(receiver). The row-normalized weight matrix W defines \(L = D - W\) and the propagator \(\exp(-\beta L)\) with β = 0.3 by default. For each target × population, the four signed shifts are summarized as a non-negative program-suppression score (total magnitude of downward shifts), which is then propagated.

Default parameters reproduce the manuscript settings: k ∈ {10, 20, 30, 50, 100}, s ∈ {0.3, 0.5, 1.0, 1.5, 2.0}, β ∈ {0.1, 0.2, 0.3, 0.5, 0.7}.

## Built-in definitions

- Four functional programs (kidney-application definitions): P1 innate sensing / NF-κB / cytokine–chemokine signaling (21 genes); P2 mitochondrial injury, oxidative stress, metabolic adaptation, adhesion, AP-1 (36 genes); P3 DNA-damage sensing, checkpoints, repair (54 genes); P4 cytosolic DNA sensing and inflammasome bridge (8 genes). Accessible as `scalepert.PROGRAMS`.
- Eighteen candidate ligand–receptor pairs for ScalePert-Tissue, including CCL5–CCR5, TNF–TNFRSF1A, HMGB1–TLR2/TLR4/AGER, VCAM1–ITGA4, VEGFA–FLT1/KDR and CXCL8–CXCR1/CXCR2. Accessible as `scalepert.LR_PAIRS`.
- Nine kidney-allograft fibrosis hub candidates as default targets: TLR2, ITGAM, CCL5, HIF1A, VCAM1, JUN, AIM2, IFI16, PYCARD. Accessible as `scalepert.DEFAULT_TARGETS`.

## Repository layout

```
src/scalepert/
├── programs.py        # program gene sets, LR pairs, gene aliases
├── preprocessing.py   # prepare_adata / score_programs
├── cell.py            # ScalePert-Cell: displacement + signed W1 shifts
├── tissue.py          # ScalePert-Tissue: LR graph + matrix-exponential propagation
├── pipeline.py        # ScalePertPipeline + parameter sensitivity
└── data.py            # example-data loader
tests/                 # pytest suite
examples/              # runnable scripts
scripts/build_example.py  # regenerates the bundled synthetic atlas
src/scalepert/data_store/ # packaged example atlas (.h5ad)
```

## Using ScalePert on your own data

Provide an `AnnData` with raw counts (or already processed values) in `.X`, a categorical population annotation in `.obs`, and your nominated target genes present in `.var_names`. Functional programs should represent biological processes relevant to your disease context; ScalePert ranks targets by how strongly their virtual reduction suppresses those programs.

```python
sens = scalepert.sensitivity(pipeline.adata_, targets=pipeline.targets_)
print(sens.ranks)          # rank of every target under each parameter setting
print(sens.correlations)   # pairwise Spearman agreement between settings
```

`sensitivity` is also importable directly: `from scalepert.pipeline import sensitivity`.

## Reproducing the manuscript analyses

The manuscript's quantitative results were produced with the research codebase used for the publication; this repository provides the reusable method implementation. To reproduce the three analytical layers:

1. **Kidney-allograft IFTA application** — run ScalePert on GSE195718 single-nucleus data (11 annotated kidney populations) with the nine default targets and the four built-in programs: see Quick start. The bundled synthetic atlas reproduces the expected qualitative behavior (IFI16 and P4 candidates prioritized) but not the published numeric values.
2. **Parameter sensitivity** — `scalepert.sensitivity` with the manuscript grids (k, s as listed above).
3. **DepMap benchmarking** — `examples/benchmark_depmap.py` with CCLE single-cell transcriptomes (≤1,000 cells per line, log1p) and DepMap Chronos gene-effect scores; reports Pearson/Spearman concordance and AUROC at Chronos < −0.5.

## Citation

If you use ScalePert, please cite the manuscript:

> ScalePert: interpretable local-manifold virtual perturbation for target prioritization from observational single-cell transcriptomes.

## License

MIT — see [LICENSE](LICENSE).
