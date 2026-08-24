# ScalePert

**ScalePert: interpretable local-manifold virtual perturbation for target prioritization from observational single-cell transcriptomes**

ScalePert is an interpretable, multiscale framework that uses naturally occurring transcriptional heterogeneity in observational single-cell (or single-nucleus) RNA-seq data to infer target-reduction-associated cellular-state transitions and prioritize candidate targets — without requiring perturbation-trained models, inferred regulatory networks, or additional molecular modalities.

The package implements the two modules described in the manuscript:

- **ScalePert-Cell** — embeds cells in a shared PCA space, builds population-restricted kNN graphs, shifts each cell toward transcriptionally similar neighbors with lower observed expression of a nominated target, reconstructs functional program scores at displaced coordinates from nearby observed cells, and quantifies the response with signed Wasserstein-1 distribution shifts.
- **ScalePert-Tissue** — maps population-level ScalePert-Cell outputs onto a ligand–receptor-derived communication graph and applies matrix-exponential (`exp(−βL)`) propagation to generate multicellular tissue-level summaries.

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
from scalepert import ScalePertPipeline, load_example_data

adata = sc.read_h5ad("your_atlas.h5ad")          # requires .obs["cell_type"]
pipeline = ScalePertPipeline(cell_type_key="cell_type", k=30, scale=1.0, beta=0.3)
pipeline.fit(adata, targets=["TLR2", "HIF1A", "AIM2", "IFI16"])

print(pipeline.cell_ranking())      # cell-population target ranking (negative = disease-program suppression)
print(pipeline.tissue_ranking())    # tissue-level ranking after graph propagation
pipeline.export("results/scalepert")
```

If your `AnnData` is already normalized, log-transformed, scaled and has `.obsm["X_pca"]` plus program scores in `.obs`, call `pipeline.fit(adata, score_existing=True)` to skip preprocessing.

### 2. Module-level usage

```python
from scalepert import prepare_adata, score_programs, scalepert_cell, build_communication_graph, propagate_tissue

adata = prepare_adata(raw_adata)                     # normalize → HVG → scale → PCA
programs = {"fibrosis": ["COL1A1", "COL3A1"], ...}  # or use built-in PROGRAMS
scored = score_programs(adata, programs)

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
import scanpy as sc
from scalepert import load_example_data, ScalePertPipeline

adata = load_example_data()   # bundled synthetic kidney-style atlas
pipe = ScalePertPipeline().fit(adata, score_existing=False)
print(pipe.cell_ranking())
```

## Method summary

Given a processed gene-by-cell matrix \(X\), population labels \(c(i)\), nominated targets and functional readouts:

1. **Shared representation** — normalize, log-transform, select highly variable genes, scale and compute PCA (50 components by default).
2. **Population-restricted neighborhoods** — cosine-distance kNN graph (k = 30 by default), constructed separately within each annotated population.
3. **Target-conditioned displacement** — for target gene *g*, the lower-expression neighbor set of cell *i* is
   \(\mathcal{N}_i^- = \{j \in \mathrm{KNN}(i) : e_j < e_i\}\),
   and the local displacement is
   \(\Delta_i = s \cdot \frac{1}{|\mathcal{N}_i^-|} \sum_{j \in \mathcal{N}_i^-} (x_j - x_i)\),
   with scaling factor \(s = 1.0\). Cells without eligible lower-expression neighbors remain at their original coordinates.
4. **Functional reconstruction** — program scores at virtual coordinates are reconstructed as the mean over nearest observed cells; cells displaced above the within-analysis 25th percentile are retained.
5. **Signed Wasserstein-1 shift** — for program *m*:
   \(\mathrm{Shift}_m = \mathrm{sign}\!\left(\bar{Q}_m - \bar{P}_m\right) \cdot W_1(P_m, Q_m)\).
   Negative values denote suppression of disease-associated programs and therefore stronger reversal priority.
6. **Additive combinations** — multi-target hypotheses average individual displacement vectors before reconstruction.
7. **Tissue aggregation** — population nodes connected by ligand–receptor-derived directed edges; row-normalized weights *W* form the propagator
   \(\exp(-\beta L),\; L = D - W,\;\beta = 0.3\) by default.

Default parameters reproduce the manuscript settings: k ∈ {10, 20, 30, 50, 100}, s ∈ {0.3, 0.5, 1.0, 1.5, 2.0}, β ∈ {0.1, 0.2, 0.3, 0.5, 0.7}.

## Repository layout

```
src/scalepert/
├── programs.py        # program gene sets, LR pairs, aliases
├── preprocessing.py   # prepare_adata / score_programs
├── cell.py            # ScalePert-Cell: displacement + signed W1 shifts
├── tissue.py          # ScalePert-Tissue: LR graph + matrix-exponential propagation
├── pipeline.py        # ScalePertPipeline + parameter sensitivity
└── data.py            # example-data loaders
tests/                 # pytest suite
examples/              # runnable scripts
data/                  # bundled example atlas
```

## Reproducing the manuscript analyses

- Kidney-allograft IFTA application: `examples/run_pipeline.py`
- Parameter sensitivity grids: `scalepert.sensitivity`
- DepMap benchmarking workflow: `examples/benchmark_depmap.py` (requires CCLE scRNA-seq and DepMap Chronos inputs)

## Citation

If you use ScalePert, please cite the manuscript:

> ScalePert: interpretable local-manifold virtual perturbation for target prioritization from observational single-cell transcriptomes.

## License

MIT — see [LICENSE](LICENSE).
