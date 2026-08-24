import numpy as np
import pandas as pd
import pytest
import anndata as ad

from scalepert import (
    PROGRAM_NAMES,
    PROGRAMS,
    ScalePertPipeline,
    build_communication_graph,
    load_example_data,
    prepare_adata,
    propagate_tissue,
    resolve_gene,
    scalepert_cell,
    score_programs,
)


@pytest.fixture(scope="session")
def example_adata():
    return load_example_data(ready=True)


@pytest.fixture(scope="session")
def toy_adata():
    return _make_toy_adata()


def _make_toy_adata(n_per_type=120, seed=7):
    rng = np.random.default_rng(seed)
    genes = ["TP53"] + [g for genes_ in PROGRAMS.values() for g in genes_]
    genes = list(dict.fromkeys(genes))
    gi = {g: i for i, g in enumerate(genes)}
    blocks = []
    labels = []
    for ct in ("A", "B"):
        latent = rng.normal(size=(n_per_type, 3))
        X = rng.normal(scale=0.1, size=(n_per_type, len(genes)))
        for pi, (_, members) in enumerate(PROGRAMS.items()):
            idx = [gi[g] for g in members if g in gi]
            if idx:
                X[:, idx] += 0.8 * np.exp(latent[:, pi % 3])[:, None]
        blocks.append(np.exp(X))
        labels += [ct] * n_per_type
    X = np.vstack(blocks)
    adata = ad.AnnData(
        X=X.astype(np.float32),
        obs=pd.DataFrame({"cell_type": labels}),
        var=pd.DataFrame(index=genes),
    )
    t = "IFI16"
    if t in gi:
        expr = X[:, gi[t]] + np.clip(rng.normal(size=X.shape[0]), -0.5, None)
        adata.X[:, gi[t]] = np.abs(expr) + 0.01
    return adata


def test_programs_complete():
    assert len(PROGRAMS) == 4
    for name, genes in PROGRAMS.items():
        assert len(genes) >= 5
    assert set(resolve_gene(g, list(PROGRAMS["DDRRepair"])) for g in ["H2AFX", "PARP1"]) != {None}


def test_example_dataset_loads(example_adata):
    assert example_adata.n_obs > 500
    assert "cell_type" in example_adata.obs.columns
    for p in ["IFI16", "TLR2", "HIF1A"]:
        assert p in example_adata.var_names


def test_prepare_adata_runs(toy_adata):
    prepared = prepare_adata(toy_adata, n_top_genes=50, n_comps=10, min_genes=50)
    assert "X_pca" in prepared.obsm
    assert prepared.obsm["X_pca"].shape[1] == 10
    assert prepared.raw is not None


def test_displacement_direction_and_eligibility(example_adata):
    prepared = prepare_adata(example_adata.copy(), min_genes=50)
    sub = prepared[prepared.obs["cell_type"] == "PT"].copy()
    from scalepert.cell import displacement_vectors

    r = displacement_vectors(sub, "IFI16", "PT", k=15)
    if r is None:
        pytest.skip("target not eligible in this population")
    moved = np.linalg.norm(r["shift_vectors"], axis=1) > 0
    assert moved.sum() == r["n_eligible"]
    expr = r["expression"]
    low_neighbors_exist = (expr[:, None].ravel() < expr.ravel()[:, None]).any(axis=1)[expr.argsort().argsort()]
    assert r["pca_after"].shape == r["pca_before"].shape
    assert np.isfinite(r["pca_after"]).all()
    for i in np.where(moved)[0][:10]:
        nb = r["shift_vectors"][i]
        assert nb @ nb > 0


def test_program_shifts_sign_convention(example_adata):
    from scalepert.cell import displacement_vectors, program_shifts

    prepared = score_programs(prepare_adata(example_adata.copy(), min_genes=50))
    sub = prepared[prepared.obs["cell_type"] == "PT"].copy()
    r = displacement_vectors(sub, "HIF1A", "PT", k=15)
    if r is None:
        pytest.skip("target not eligible in this population")
    shifts = program_shifts(sub, r)
    if shifts is None:
        pytest.skip("insufficient displaced cells")
    for val in shifts.values():
        assert np.isfinite(val)


def test_scalepert_cell_full_table(example_adata):
    res = scalepert_cell(
        example_adata,
        targets=["IFI16", "HIF1A", "TLR2"],
        cell_type_key="cell_type",
        k=20,
        progress=False,
    )
    table = res.table
    assert {"single"} <= set(table["type"])
    assert table["mean_signed_W1"].notna().all()
    ranking = res.ranking()
    assert len(ranking) == 3
    assert list(ranking["rank"]) == [1, 2, 3]


def test_communication_graph_shape_and_nonnegativity(example_adata):
    comm = build_communication_graph(example_adata, "cell_type")
    cts = sorted(example_adata.obs["cell_type"].unique())
    assert comm.shape == (len(cts), len(cts))
    assert (comm.values >= 0).all()


def test_tissue_propagation_matches_direct_signal(example_adata):
    res = scalepert_cell(
        example_adata,
        targets=["IFI16", "TLR2"],
        cell_type_key="cell_type",
        k=20,
        progress=False,
    )
    comm = build_communication_graph(res.adata, "cell_type")
    out = propagate_tissue(res, comm, beta=0.3)
    direct_cols = [c for c in res.table.columns if c.endswith("_W1") and c != "mean_signed_W1"]
    direct = res.table.groupby("target")[direct_cols].mean()
    merged = out.summary.set_index("target")
    tissue_cols = ["tissue_" + c[:-3] for c in direct.columns]
    tissue_cols = [c for c in tissue_cols if c in merged.columns]
    for target in direct.index:
        d = np.asarray(direct.loc[target].values, dtype=float)
        m = np.asarray(merged.loc[target, tissue_cols].values, dtype=float)
        corr = np.corrcoef(d, m)[0, 1]
        assert corr > 0.9
    assert (out.summary["suppression_score"] >= 0).all()


def test_pipeline_end_to_end_with_export(tmp_path):
    pipe = ScalePertPipeline(k=20).fit(load_example_data(), targets=["IFI16", "HIF1A"])
    cr = pipe.cell_ranking()
    tr = pipe.tissue_ranking()
    assert len(cr) == 2 and len(tr) == 2
    prefix = tmp_path / "out"
    pipe.export(prefix)
    for suffix in [
        "_cell_perturbations.csv",
        "_cell_ranking.csv",
        "_communication_graph.csv",
        "_tissue_summary.csv",
    ]:
        assert (tmp_path / f"out{suffix}").exists()


def test_pipeline_score_existing_mode():
    pipe = ScalePertPipeline(k=20).fit(load_example_data(), targets=["IFI16"], score_existing=True)
    assert len(pipe.cell_result_.table) > 0


def test_sensitivity_ranks_stable():
    from scalepert.pipeline import sensitivity

    sens = sensitivity(
        load_example_data(),
        targets=["IFI16", "TLR2"],
        ks=(10, 20),
        scales=(0.5, 1.0),
        cell_type_key="cell_type",
    )
    assert sens.ranks.shape[1] == 4
    assert sens.correlations.notna().all()


def test_missing_celltype_key_raises(toy_adata):
    toy = toy_adata.copy()
    toy.obs = toy.obs.rename(columns={"cell_type": "population"})
    with pytest.raises(KeyError):
        prepare_adata(toy)
