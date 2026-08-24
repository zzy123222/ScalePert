import numpy as np
import pandas as pd
import anndata as ad

from scalepert.programs import resolve_programs, DEFAULT_TARGETS, HUB_PROGRAM, PROGRAM_NAMES
from scalepert.data import _default_gene_pool


TARGET_DISEASE_COUPLING = {
    "IFI16": 1.00,
    "AIM2": 0.85,
    "PYCARD": 0.80,
    "HIF1A": 0.55,
    "VCAM1": 0.50,
    "JUN": 0.45,
    "TLR2": 0.35,
    "CCL5": 0.30,
    "ITGAM": 0.25,
}

CELL_TYPE_WEIGHTS = {"IMM": 1.15, "FIB": 1.05, "PT": 0.90}


def _make_latent(rng, n):
    latent = rng.normal(size=(n, PROGRAM_NAMES.__len__()))
    return latent


def _program_expression(rng, severity, gene_index, programs):
    n = len(severity)
    base = rng.normal(scale=0.08, size=(n, len(gene_index)))
    for pi, (_, genes) in enumerate(programs.items()):
        members = [gene_index[g] for g in genes if g in gene_index]
        if not members:
            continue
        base[:, members] += (0.9 + 0.25 * pi) * severity[:, None]
    return base


def _target_modulation(severity, gene_index):
    mod = {}
    for t, strength in TARGET_DISEASE_COUPLING.items():
        if t not in gene_index:
            continue
        mod[t] = strength * severity + np.random.default_rng(7).normal(
            scale=0.12, size=len(severity)
        )
    return mod


def build_example_atlas(
    n_per_type=400,
    cell_types=("PT", "FIB", "IMM"),
    seed=42,
):
    rng = np.random.default_rng(seed)
    genes = _default_gene_pool()
    gene_index = {g: i for i, g in enumerate(genes)}
    blocks = []
    obs_rows = []
    for ci, ct in enumerate(cell_types):
        weight = CELL_TYPE_WEIGHTS.get(ct, 1.0)
        disease_state = "IFTA" if ct != "PT" else ("IFTA" if ci % 2 == 0 else "stable")
        frac_diseased = 0.75 if disease_state == "IFTA" else 0.25
        diseased = rng.random(n_per_type) < frac_diseased
        severity = np.clip(0.4 * weight * rng.normal(size=n_per_type) + 1.1 * weight * diseased, -1.5, None)
        base = _program_expression(rng, severity, gene_index, resolve_programs(genes))
        counts = np.exp(base)
        mod = _target_modulation(severity, gene_index)
        for t, m in mod.items():
            idx = gene_index.get(t)
            if idx is not None:
                counts[:, idx] = np.abs(counts[:, idx] + m) + 1e-3
        blocks.append(counts)
        obs_rows.append(
            pd.DataFrame(
                {
                    "cell_type": [ct] * n_per_type,
                    "disease": ["IFTA" if d else "stable" for d in diseased],
                    "batch": [f"donor_{ci}"] * n_per_type,
                    "severity": severity.astype(np.float32),
                }
            )
        )

    X = np.vstack(blocks).astype(np.float32)
    obs = pd.concat(obs_rows, ignore_index=True)
    var = pd.DataFrame(index=pd.Index(genes))
    adata = ad.AnnData(X=X, obs=obs, var=var)
    adata.obs_names = [f"cell_{i:05d}" for i in range(adata.n_obs)]
    adata.obs["cell_type"] = adata.obs["cell_type"].astype("category")
    return adata


def main(path=None):
    adata = build_example_atlas()
    from scalepert.preprocessing import prepare_adata, score_programs

    prepared = prepare_adata(adata, min_genes=50)
    prepared = score_programs(prepared)
    if path is None:
        import os

        path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "src", "scalepert", "data_store", "scalepert_example.h5ad")
    prepared.write_h5ad(path)
    return path


if __name__ == "__main__":
    print(main())
