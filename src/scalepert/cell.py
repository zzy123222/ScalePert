import numpy as np
import pandas as pd
from scipy.stats import wasserstein_distance
from sklearn.neighbors import NearestNeighbors

from .programs import resolve_gene, resolve_programs


def displacement_vectors(
    adata,
    target,
    cell_type,
    cell_type_key="cell_type",
    k=30,
    scale=1.0,
    min_cells=50,
    pca_key="X_pca",
):
    mask = (adata.obs[cell_type_key] == cell_type).values
    idx = np.where(mask)[0]
    if len(idx) < min_cells:
        return None
    gene = resolve_gene(target, adata.var_names)
    if gene is None:
        return None
    X_pca = np.asarray(adata.obsm[pca_key][idx])
    g_idx = list(adata.var_names).index(gene)
    raw = adata.X[idx]
    if hasattr(raw, "toarray"):
        raw = raw.toarray()
    expr = raw[:, g_idx].astype(float)

    n_neighbors = min(k, len(idx) - 1)
    if n_neighbors < 1:
        return None
    nn = NearestNeighbors(n_neighbors=n_neighbors, metric="cosine")
    nn.fit(X_pca)
    _, nbrs = nn.kneighbors(X_pca)

    lower = expr[nbrs] < expr[:, None]
    counts = lower.sum(axis=1)
    eligible = counts > 0
    if not eligible.any():
        return None

    sums = np.zeros((len(idx), X_pca.shape[1]))
    for i in range(len(idx)):
        nb = nbrs[i][lower[i]]
        if len(nb) > 0:
            sums[i] = X_pca[nb].mean(axis=0) - X_pca[i]
    shift = sums * scale
    pca_after = X_pca + shift
    return {
        "target": target,
        "gene": gene,
        "cell_type": cell_type,
        "cell_indices": idx,
        "pca_before": X_pca,
        "pca_after": pca_after,
        "shift_vectors": shift,
        "expression": expr,
        "n_cells": int(len(idx)),
        "n_eligible": int(eligible.sum()),
        "obs": adata.obs.iloc[idx],
    }


def program_shifts(
    adata,
    result,
    programs=None,
    reconstruction_neighbors=10,
    moved_percentile=25,
    min_moved=10,
):
    if result is None:
        return None
    if programs is None:
        programs = resolve_programs(adata.var_names)
    obs = result["obs"]
    scores = {m: obs[m].values.astype(float) for m in programs if m in obs.columns}
    if not scores:
        return None
    mag = np.linalg.norm(result["shift_vectors"], axis=1)
    threshold = np.percentile(mag, moved_percentile)
    moved = mag > threshold
    if moved.sum() < min_moved:
        return None
    pca_before = result["pca_before"]
    pca_after = result["pca_after"]
    n_post = min(reconstruction_neighbors, len(scores[next(iter(scores))]) - 1)
    if n_post < 1:
        return None
    nn = NearestNeighbors(n_neighbors=n_post, metric="euclidean")
    nn.fit(pca_before)
    _, nbrs_post = nn.kneighbors(pca_after)
    shifts = {}
    for m, s in scores.items():
        after = s[nbrs_post].mean(axis=1)
        w1 = wasserstein_distance(s[moved], after[moved])
        direction = np.sign(after[moved].mean() - s[moved].mean())
        shifts[m] = float(direction * w1)
    result = dict(shifts)
    result["n_cells"] = result.get("n_cells", result.get("n_cells"))
    return shifts


def scalepert_cell(
    adata,
    targets,
    cell_types=None,
    cell_type_key="cell_type",
    k=30,
    scale=1.0,
    min_cells=50,
    programs=None,
    combos=None,
    pca_key="X_pca",
    progress=True,
):
    from .programs import DEFAULT_TARGETS, HUB_PROGRAM, resolve_programs

    targets = list(targets) if targets is not None else list(DEFAULT_TARGETS)
    if cell_types is None:
        cell_types = sorted(adata.obs[cell_type_key].astype(str).unique())
    if programs is None:
        obs_programs = [m for m in resolve_programs(adata.var_names) if m in adata.obs.columns]
        programs = {m: [] for m in obs_programs} if obs_programs else resolve_programs(adata.var_names)
    rows = []
    results_store = {}

    def run_one(gene_set, label, kind, weight=None):
        for ct in cell_types:
            per_gene_shifts = []
            eligible_flags = []
            n_cells_ct = None
            for g in gene_set:
                r = displacement_vectors(
                    adata,
                    g,
                    ct,
                    cell_type_key=cell_type_key,
                    k=k,
                    scale=weight if weight is not None else scale,
                    min_cells=min_cells,
                    pca_key=pca_key,
                )
                if r is not None:
                    n_cells_ct = r["n_cells"]
                sh = program_shifts(adata, r, programs=programs)
                if sh is not None and all(np.isfinite(v) for v in sh.values()):
                    per_gene_shifts.append(sh)
                    eligible_flags.append(True)
                else:
                    eligible_flags.append(False)
            if not per_gene_shifts:
                continue
            keys = list(per_gene_shifts[0].keys())
            combined = {kk: float(np.mean([sh[kk] for sh in per_gene_shifts])) for kk in keys}
            row = {
                "target": label,
                "type": kind,
                "cell_type": ct,
                "n_genes": len(gene_set),
                "program": "+".join(sorted({HUB_PROGRAM.get(g, "?") for g in gene_set})),
                "n_cells": n_cells_ct,
                "eligible_targets": sum(eligible_flags),
            }
            for m in programs:
                row[f"{m}_W1"] = combined.get(m, np.nan)
            row["mean_signed_W1"] = float(np.mean(list(combined.values())))
            neg = [v for v in combined.values() if v < 0]
            pos = [v for v in combined.values() if v > 0]
            row["reversal_score"] = float(sum(-v for v in neg))
            row["upshift_score"] = float(sum(pos))
            rows.append(row)

    for t in targets:
        if progress:
            print(f"[ScalePert-Cell] {t}")
        run_one([t], t, "single")
    if combos:
        for combo in combos:
            label = "+".join(combo)
            if progress:
                print(f"[ScalePert-Cell] combo {label}")
            run_one(list(combo), label, "combo", weight=scale / max(len(combo), 1))

    df = pd.DataFrame(rows)
    if len(df) == 0:
        raise RuntimeError(
            "no eligible perturbations produced; check that targets are expressed "
            "and populations contain at least min_cells cells"
        )
    single = df[df["type"] == "single"].copy()
    if len(single) > 0:
        rank_values = single.groupby("target")["mean_signed_W1"].mean().sort_values()
        df["rank"] = df["target"].map({t: i + 1 for i, t in enumerate(rank_values.index)})
    else:
        df["rank"] = np.arange(1, len(df) + 1)
    return CellResult(df, adata)


class CellResult:
    def __init__(self, table, adata):
        self.table = table
        self.adata = adata

    def __repr__(self):
        return f"CellResult(perturbations={len(self.table)})"

    def ranking(self, score="mean_signed_W1", ascending=True):
        single = self.table[self.table["type"] == "single"]
        agg = (
            single.groupby("target")[score]
            .mean()
            .sort_values(ascending=ascending)
            .reset_index()
        )
        agg.insert(0, "rank", np.arange(1, len(agg) + 1))
        return agg

    def tissue_summary(self, score="mean_signed_W1", ascending=True):
        sub = self.table[self.table["type"] == "single"]
        out = (
            sub.groupby("target")[score]
            .agg(["mean", "median", "min", "max", "count"])
            .sort_values("mean", ascending=ascending)
            .reset_index()
        )
        out.columns = ["target", "mean_score", "median_score", "min_score", "max_score", "n_populations"]
        return out
