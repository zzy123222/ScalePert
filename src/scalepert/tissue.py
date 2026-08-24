import numpy as np
import pandas as pd
from scipy.linalg import expm

from .programs import LR_PAIRS, GENE_ALIASES


def _resolve_lr_pairs(var_names):
    names = set(map(str, var_names))

    def present(gene):
        if gene in names:
            return gene
        alias = GENE_ALIASES.get(gene)
        if alias and alias in names:
            return alias
        return None

    pairs = []
    for lig, rec in LR_PAIRS:
        l = present(lig)
        r = present(rec)
        if l and r:
            pairs.append((l, r))
    return pairs


def build_communication_graph(adata, cell_type_key="cell_type", lr_pairs=None):
    if lr_pairs is not None:
        var_names = set(map(str, adata.var_names))
        pairs = [(l, r) for l, r in lr_pairs if l in var_names and r in var_names]
    else:
        pairs = _resolve_lr_pairs(adata.var_names)
    if not pairs:
        raise ValueError(
            "none of the ligand-receptor pairs are present in the matrix; "
            "supply custom lr_pairs matching your gene names"
        )
    X = adata.X
    if hasattr(X, "toarray"):
        X = X.toarray()
    X = np.asarray(X)
    cell_types = sorted(adata.obs[cell_type_key].astype(str).unique())
    g2i = {g: i for i, g in enumerate(adata.var_names)}
    labels = adata.obs[cell_type_key].astype(str).values
    expr = {g: X[:, g2i[g]].astype(float) for pair in pairs for g in pair}
    W = np.zeros((len(cell_types), len(cell_types)))
    for i, sender in enumerate(cell_types):
        ms = labels == sender
        for j, receiver in enumerate(cell_types):
            mr = labels == receiver
            total = 0.0
            for lig, rec in pairs:
                total += expr[lig][ms].mean() * expr[rec][mr].mean()
            W[i, j] = total
    return pd.DataFrame(W, index=cell_types, columns=cell_types)


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
    detail_rows = []
    for label, sub in table.groupby("target"):
        ptype = sub["type"].iloc[0]
        vec_by_program = {}
        for p in programs:
            col = f"{p}_W1"
            v = []
            for ct in cts:
                hit = sub[sub["cell_type"] == ct]
                v.append(float(hit[col].values[0]) if len(hit) else 0.0)
            vec_by_program[p] = np.asarray(v)
        suppression_per_ct = np.zeros(len(cts))
        prop_means = {}
        prop_by_program = {}
        for p in programs:
            prop = operator @ vec_by_program[p]
            prop_by_program[p] = prop
            prop_means[p] = float(prop.mean())
            suppression_per_ct += np.maximum(-prop, 0.0)
        row = {
            "target": label,
            "type": ptype,
            "suppression_score": float(suppression_per_ct.sum()),
            "propagated_signed_mean": float(np.mean(list(prop_means.values()))),
        }
        for p in programs:
            row[f"tissue_{p}"] = prop_means[p]
        rows.append(row)
        for i, ct in enumerate(cts):
            d = {"target": label, "cell_type": ct}
            for p in programs:
                d[f"{p}_W1"] = float(prop_by_program[p][i])
            d["suppression_score"] = float(suppression_per_ct[i])
            detail_rows.append(d)
    return TissueResult(pd.DataFrame(rows), pd.DataFrame(detail_rows), beta)


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
