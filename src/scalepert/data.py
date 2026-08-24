import os

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
