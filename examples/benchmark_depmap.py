import argparse
import numpy as np
import pandas as pd
import scanpy as sc
from scipy.stats import spearmanr, pearsonr
from sklearn.metrics import roc_auc_score

from scalepert import prepare_adata, score_programs, scalepert_cell


def load_inputs(adata_path, chronos_path):
    adata = sc.read_h5ad(adata_path)
    chronos = pd.read_csv(chronos_path, index_col=0)
    return adata, chronos


def run_scalepert_ranking(adata, targets, cell_type_key="cell_type"):
    prepared = prepare_adata(adata)
    scored = score_programs(prepared)
    res = scalepert_cell(scored, targets=targets, cell_type_key=cell_type_key)
    ranking = res.ranking()
    return ranking


def benchmark(ranking, chronos, essentiality_threshold=-0.5):
    merged = ranking.merge(
        chronos.rename(columns={chronos.columns[0]: "chronos"}), left_on="target", right_index=True
    )
    essential = (merged["chronos"] < essentiality_threshold).astype(int)
    pearson_r, pearson_p = pearsonr(merged["mean_signed_W1"], merged["chronos"])
    spearman_r, spearman_p = spearmanr(merged["mean_signed_W1"], merged["chronos"])
    auc = roc_auc_score(essential, -merged["mean_signed_W1"])
    return {
        "n_genes": len(merged),
        "pearson_r": pearson_r,
        "pearson_p": pearson_p,
        "spearman_rho": spearman_r,
        "spearman_p": spearman_p,
        "pooled_auroc": auc,
    }


def main():
    parser = argparse.ArgumentParser(description="ScalePert vs DepMap Chronos benchmark")
    parser.add_argument("adata", help="h5ad atlas of cell lines (one population per line)")
    parser.add_argument("chronos", help="CSV with gene x line Chronos scores")
    parser.add_argument("--targets", nargs="+", required=True)
    parser.add_argument("--cell-type-key", default="cell_type")
    args = parser.parse_args()

    adata, chronos = load_inputs(args.adata, args.chronos)
    ranking = run_scalepert_ranking(adata, args.targets, args.cell_type_key)
    metrics = benchmark(ranking, chronos)
    print(pd.Series(metrics).to_string())


if __name__ == "__main__":
    main()
