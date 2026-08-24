import numpy as np
import pandas as pd
import anndata as ad
from scipy.stats import wasserstein_distance

from scalepert.programs import resolve_programs, DEFAULT_TARGETS
from scalepert.data import _default_gene_pool


def _make_latent(rng, n):
    return rng.normal(size=(n, 4))


def _program_expression(rng, latent, gene_index, programs):
    n, p = len(latent), len(gene_index)
    base = rng.normal(scale=0.05, size=(n, p))
    for pi, (_, genes) in enumerate(programs.items()):
        members = [gene_index[g] for g in genes if g in gene_index]
        if not members:
            continue
        factor = latent[:, pi % latent.shape[1]]
        boost = 0.6 * np.exp(0.8 * factor)[:, None]
        base[:, members] += boost
    return base


def _target_modulation(rng, latent, gene_index, targets):
    mod = {}
    for ti, t in enumerate(targets):
        if t not in gene_index:
            continue
        direction = 1.0 if ti % 2 == 0 else -1.0
        mod[t] = direction * (0.5 + 0.2 * latent[:, ti % latent.shape[1]])
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
        latent = _make_latent(rng, n_per_type)
        base = _program_expression(rng, latent, gene_index, resolve_programs(genes))
        counts = np.exp(base)
        mod = _target_modulation(rng, latent, gene_index, DEFAULT_TARGETS)
        blocks.append(counts)
        obs_rows.append(
            pd.DataFrame(
                {
                    "cell_type": [ct] * n_per_type,
                    "disease": ["IFTA" if ci % 2 == 0 else "stable"] * n_per_type,
                    "batch": [f"donor_{ci}"] * n_per_type,
                }
            )
        )
        for t, m in mod.items():
            idx = gene_index.get(t)
            if idx is not None:
                counts[:, idx] = np.abs(counts[:, idx] + np.clip(m, -1.5, None))

    X = np.vstack(blocks).astype(np.float32)
    obs = pd.concat(obs_rows, ignore_index=True)
    var = pd.DataFrame(index=pd.Index(genes))
    adata = ad.AnnData(X=X, obs=obs, var=var)
    adata.obs_names = [f"cell_{i:05d}" for i in range(adata.n_obs)]
    adata.obs["cell_type"] = adata.obs["cell_type"].astype("category")

    for name, genes_ in resolve_programs(adata.var_names).items():
        members = [g for g in genes_ if g in set(adata.var_names)]
        scores = np.zeros(adata.n_obs)
        if len(members) >= 3:
            sub = adata[:, members].X
            sub = np.asarray(sub) if not hasattr(sub, "toarray") else sub.toarray()
            ref = adata.X
            ref = np.asarray(ref) if not hasattr(ref, "toarray") else ref.toarray()
            scores = sub.mean(axis=1) - ref.mean(axis=1)
        adata.obs[name] = scores.astype(float)

    return adata


def write_example_dataset(path=None):
    adata = build_example_atlas()
    if path is None:
        path = _default_store()
    adata.write_h5ad(path)
    return path


def _default_store():
    import os

    return os.path.join(os.path.dirname(os.path.abspath(__file__)), "data_store", "scalepert_example.h5ad")


if __name__ == "__main__":
    print(write_example_dataset())
