import numpy as np
import scanpy as sc


def prepare_adata(
    adata,
    n_top_genes=3000,
    n_comps=50,
    n_neighbors=30,
    cell_type_key="cell_type",
    min_genes=200,
    store_normalized=True,
):
    if cell_type_key not in adata.obs.columns:
        raise KeyError(
            f"cell type annotation column '{cell_type_key}' not found; "
            "pass cell_type_key explicitly"
        )
    adata = adata.copy()
    detected = _genes_per_cell(adata)
    if adata.n_vars >= min_genes and np.median(detected) >= min_genes:
        sc.pp.filter_cells(adata, min_genes=min_genes)
    sc.pp.filter_genes(adata, min_cells=3)
    if adata.n_obs == 0 or adata.n_vars == 0:
        raise ValueError(
            "filtering removed all cells/genes; lower min_genes/min_cells or check counts"
        )
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    if store_normalized:
        adata.raw = adata
    max_val = adata.X.max() if hasattr(adata.X, "max") else 0
    flavor = "seurat_v3" if max_val > 20 else "seurat"
    work = adata.copy()
    sc.pp.highly_variable_genes(work, n_top_genes=min(n_top_genes, work.n_vars), flavor=flavor)
    hvg_idx = work.var_names[work.var["highly_variable"]]
    sub = adata[:, hvg_idx].copy()
    sc.pp.scale(sub, max_value=10)
    n_comp = min(n_comps, sub.n_vars - 1, sub.n_obs - 1)
    sc.tl.pca(sub, n_comps=n_comp)
    adata.obsm["X_pca"] = sub.obsm["X_pca"]
    sc.pp.neighbors(adata, n_neighbors=min(n_neighbors, adata.n_obs - 1), use_rep="X_pca")
    return adata


def _genes_per_cell(adata):
    X = adata.X
    if hasattr(X, "toarray"):
        X = X.toarray()
    return np.asarray((np.asarray(X) > 0).sum(axis=1)).ravel()


def score_programs(adata, programs=None, ctrl_as_ref=False):
    from .programs import resolve_programs

    if programs is None:
        programs = resolve_programs(adata.var_names)
    adata = adata.copy()
    for name, genes in programs.items():
        present = [g for g in genes if g in set(adata.var_names)]
        if len(present) >= 3:
            _score_one(adata, name, present, ctrl_as_ref)
        else:
            adata.obs[name] = 0.0
    return adata


def _score_one(adata, name, present, ctrl_as_ref):
    try:
        work = adata.raw.to_adata() if adata.raw is not None and adata.raw.n_vars > adata.n_vars else adata
        sc.tl.score_genes(
            work,
            gene_list=present,
            score_name=name,
            use_raw=False,
            ctrl_as_ref=ctrl_as_ref,
        )
        adata.obs[name] = work.obs[name].values
    except (RuntimeError, ValueError):
        X = adata[:, present].X
        if hasattr(X, "toarray"):
            X = X.toarray()
        X = np.asarray(X, dtype=float)
        z = (X - X.mean(axis=0)) / (X.std(axis=0) + 1e-9)
        adata.obs[name] = z.mean(axis=1)
