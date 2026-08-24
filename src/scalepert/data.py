import os

import numpy as np
import pandas as pd
import anndata as ad


def load_example_data(ready=True):
    path = _example_path("scalepert_example.h5ad")
    adata = ad.read_h5ad(path)
    if not ready:
        return adata
    from .preprocessing import prepare_adata, score_programs
    from .programs import resolve_programs

    prepared = prepare_adata(adata, min_genes=50)
    prepared = score_programs(prepared, resolve_programs(prepared.var_names))
    return prepared


def _example_path(name):
    return os.path.join(os.path.dirname(os.path.abspath(__file__)), "data_store", name)


def simulate_atlas(
    n_per_type=400,
    cell_types=("Alpha", "Beta", "Ductal"),
    genes=None,
    programs=None,
    targets=None,
    seed=42,
):
    rng = np.random.default_rng(seed)
    genes = list(genes) if genes is not None else _default_gene_pool()
    programs = programs or PROGRAMS
    targets = list(targets) if targets is not None else DEFAULT_TARGETS

    gene_index = {g: i for i, g in enumerate(genes)}
    blocks = []
    obs_rows = []
    for ci, ct in enumerate(cell_types):
        latent = rng.normal(size=(n_per_type, 4))
        base = rng.normal(scale=0.05, size=(n_per_type, len(genes)))
        for pi, (pname, pgenes) in enumerate(programs.items()):
            members = [gene_index[g] for g in pgenes if g in gene_index]
            if not members:
                continue
            base[:, members] += 0.6 * latent[:, pi % latent.shape[1]][:, None]
        block = np.exp(base)
        blocks.append(block)
        obs_rows.append(pd.DataFrame({"cell_type": ct, "batch": f"b{ci}", "index": np.arange(n_per_type)}))

    X = np.vstack(blocks)
    var = pd.DataFrame(index=genes)
    obs = pd.concat(obs_rows, ignore_index=True)
    adata = ad.AnnData(X=X, obs=obs, var=var)
    adata.obs_names = [f"cell_{i}" for i in range(adata.n_obs)]
    adata.var_names = pd.Index(genes)

    target_expr = {}
    for t in targets:
        if t in gene_index:
            col = gene_index[t]
            modulation = 0.8 * latent_modulation(adata.n_obs, rng)
            vals = X[:, col]
            vals = np.clip(vals + modulation, 0, None)
            X[:, col] = vals
            target_expr[t] = True
    adata.X = X
    return adata


def latent_modulation(n, rng):
    return rng.normal(loc=0.0, scale=0.15, size=n)


def _default_gene_pool():
    from .programs import PROGRAMS

    pool = []
    for pgenes in PROGRAMS.values():
        pool.extend(pgenes)
    seen = set()
    out = []
    for g in pool:
        if g not in seen:
            seen.add(g)
            out.append(g)
    extras = ["GAPDH", "ACTB", "B2M", "TMSB4X", "FTL", "FTH1"]
    for g in extras:
        if g not in seen:
            seen.add(g)
            out.append(g)
    return out
