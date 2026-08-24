import numpy as np
import pandas as pd

from .cell import scalepert_cell
from .preprocessing import prepare_adata, score_programs
from .programs import DEFAULT_TARGETS, resolve_programs
from .tissue import build_communication_graph, propagate_tissue


def sensitivity(
    adata,
    targets=None,
    ks=(10, 20, 30, 50, 100),
    scales=(0.3, 0.5, 1.0, 1.5, 2.0),
    cell_type_key="cell_type",
    programs=None,
):
    from scipy.stats import spearmanr

    if targets is None:
        targets = DEFAULT_TARGETS
    rankings = {}
    for k in ks:
        res = scalepert_cell(
            adata,
            targets,
            cell_type_key=cell_type_key,
            k=k,
            programs=programs,
            progress=False,
        )
        rankings[f"k={k}"] = res.ranking().set_index("target")["rank"]
    for s in scales:
        res = scalepert_cell(
            adata,
            targets,
            cell_type_key=cell_type_key,
            k=30,
            scale=s,
            programs=programs,
            progress=False,
        )
        rankings[f"s={s}"] = res.ranking().set_index("target")["rank"]
    grid = pd.DataFrame(rankings)
    corr = {}
    keys = list(grid.columns)
    for i, a in enumerate(keys):
        for b in keys[i + 1 :]:
            rho, _ = spearmanr(grid[a], grid[b])
            corr[f"{a}|{b}"] = rho
    return SensitivityResult(grid, pd.Series(corr, name="spearman"))


class SensitivityResult:
    def __init__(self, ranks, correlations):
        self.ranks = ranks
        self.correlations = correlations

    def __repr__(self):
        return f"SensitivityResult(parameters={len(self.ranks.columns)})"


class ScalePertPipeline:
    def __init__(
        self,
        cell_type_key="cell_type",
        k=30,
        scale=1.0,
        min_cells=50,
        beta=0.3,
        n_top_genes=3000,
        n_comps=50,
        lr_pairs=None,
    ):
        self.cell_type_key = cell_type_key
        self.k = k
        self.scale = scale
        self.min_cells = min_cells
        self.beta = beta
        self.n_top_genes = n_top_genes
        self.n_comps = n_comps
        self.lr_pairs = lr_pairs

    def fit(self, adata, targets=None, combos=None, score_existing=False):
        self.targets_ = list(targets) if targets is not None else list(DEFAULT_TARGETS)
        if score_existing:
            work = adata.copy()
            self.programs_ = resolve_programs(work.var_names)
            work = score_programs(work, self.programs_)
        else:
            work = prepare_adata(adata, n_top_genes=self.n_top_genes, n_comps=self.n_comps)
            self.programs_ = resolve_programs(work.var_names)
            work = score_programs(work, self.programs_)
        self.adata_ = work
        self.cell_result_ = scalepert_cell(
            work,
            self.targets_,
            cell_type_key=self.cell_type_key,
            k=self.k,
            scale=self.scale,
            min_cells=self.min_cells,
            programs=self.programs_,
            combos=combos,
        )
        try:
            self.comm_ = build_communication_graph(work, self.cell_type_key, lr_pairs=self.lr_pairs)
            self.tissue_result_ = propagate_tissue(self.cell_result_, self.comm_, beta=self.beta)
        except ValueError:
            self.comm_ = None
            self.tissue_result_ = None
        return self

    def cell_ranking(self):
        return self.cell_result_.ranking()

    def tissue_ranking(self):
        if self.tissue_result_ is None:
            raise RuntimeError("communication graph unavailable; tissue ranking not computed")
        return self.tissue_result_.ranking()

    def summary_table(self):
        rows = [self.cell_result_.tissue_summary()]
        rows[0].insert(0, "layer", "ScalePert-Cell")
        out = rows[0]
        if self.tissue_result_ is not None:
            t = self.tissue_result_.summary[self.tissue_result_.summary["type"] == "single"][
                ["target", "suppression_score"]
            ].copy()
            t.insert(0, "layer", "ScalePert-Tissue")
            out = pd.concat([out, t], ignore_index=True)
        return out

    def export(self, path_prefix):
        base = str(path_prefix)
        self.cell_result_.table.to_csv(f"{base}_cell_perturbations.csv", index=False)
        self.cell_ranking().to_csv(f"{base}_cell_ranking.csv", index=False)
        if self.comm_ is not None:
            self.comm_.to_csv(f"{base}_communication_graph.csv")
        if self.tissue_result_ is not None:
            self.tissue_result_.summary.to_csv(f"{base}_tissue_summary.csv", index=False)
            self.tissue_result_.per_population.to_csv(
                f"{base}_tissue_per_population.csv", index=False
            )
