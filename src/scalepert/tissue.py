import numpy as np
import pandas as pd
from scipy.linalg import expm

from .programs import LR_PAIRS


def build_communication_graph(adata, cell_type_key="cell_type", lr_pairs=None, use_layer=None):
    if lr_pairs is None:
        lr_pairs = LR_PAIRS
    var_names = set(map(str, adata.var_names))
    pairs = [(l, r) for l, r in lr_pairs if l in var_names and r in var_names]
    if not pairs:
        raise ValueError(
            "none of the ligand-receptor pairs are present in the matrix; "
            "supply custom lr_pairs matching your gene names"
        )
    X = adata.X
    if hasattr(X, "toarray"):
        X = X.toarray()
    cell_types = sorted(adata.obs[cell_type_key].astype(str).unique())
    g2i = {g: i for i, g in enumerate(adata.var_names)}
    labels = adata.obs[cell_type_key].astype(str).values
    expr = {g: X[:, g2i[g]].astype(float) for l, r in pairs for g in (l, r)}
    W = np.zeros((len(cell_types), len(cell_types)))
    for i, sender in enumerate(cell_types):
        ms = labels == sender
        for j, receiver in enumerate(cell_types):
            mr = labels == receiver
            total = 0.0
            for lig, rec in pairs:
                total += expr[lig][ms].mean() * expr[rec][mr].mean()
            W[i, j] = total
    comm = pd.DataFrame(W, index=cell_types, columns=cell_types)
    return comm


def propagate_tissue(cell_result, comm, beta=0.3):
    table = cell_result.table if hasattr(cell_result, "table") else cell_result
    program_cols = [c for c in table.columns if c.endswith("_W1") and c != "mean_signed_W1"]
    programs = [c[:-3] for c in program_cols]
    cts = list(comm.index)
    W = comm.values.astype(float).copy()
    row_sums = W.sum(axis=1)
    row_sums[row_sums == 0] = 1.0
    W = W / row_sums[:, None]
    L = np.diag(W.sum(axis=1)) - W
    operator = expm(-beta * L)
    rows = []
    for label, sub in table.groupby("target"):
        ptype = sub["type"].iloc[0]
        vec = {p: [] for p in programs}
        for ct in cts:
            hit = sub[sub["cell_type"] == ct]
            for p in programs:
                col = f"{p}_W1"
                vec[p].append(float(hit[col].values[0]) if len(hit) and col in hit.columns else 0.0)
        row = {"target": label, "type": ptype}
        propagated_means = {}
        propagated_by_ct = {}
        for p in programs:
            v = np.asarray(vec[p])
            prop = operator @ v
            propagated_means[p] = float(prop.mean())
            for i, ct in enumerate(cts):
                propagated_by_ct[(ct, p)] = float(prop[i])
        for p in programs:
            row[f"tissue_{p}"] = propagated_means[p]
        neg = [v for v in propagated_means.values() if v < 0]
        row["suppression_score"] = float(sum(-v for v in neg))
        row["propagated_signed_mean"] = float(np.mean(list(propagated_means.values())))
        rows.append(row)
    out = pd.DataFrame(rows)
    detail_rows = []
    for _, r in out.iterrows():
        for ct in cts:
            d = {"target": r["target"], "cell_type": ct}
            for p in programs:
                d[f"{p}_W1"] = propagated_by_ct.get((ct, p), np.nan)
            detail_rows.append(d)
    return TissueResult(out, pd.DataFrame(detail_rows), beta)


class TissueResult:
    def __init__(self, summary, per_population, beta):
        self.summary = summary
        self.per_population = per_population
        self.beta = beta

    def __repr__(self):
        return f"TissueResult(targets={len(self.summary)}, beta={self.beta})"

    def ranking(self, score="suppression_score", ascending=False):
        single = self.summary[self.summary["type"] == "single"]
        return (
            single.sort_values(score, ascending=ascending)[["target", score]]
            .reset_index(drop=True)
            .rename(columns={score: "score"})
        )
