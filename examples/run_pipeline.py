import scanpy as sc
from scalepert import load_example_data, ScalePertPipeline


def main():
    adata = load_example_data()
    print(f"example atlas: {adata.shape[0]} cells x {adata.shape[1]} genes")
    print(adata.obs["cell_type"].value_counts())

    pipeline = ScalePertPipeline(cell_type_key="cell_type", k=30, beta=0.3)
    pipeline.fit(adata)

    print("\nScalePert-Cell ranking (negative = stronger program suppression):")
    print(pipeline.cell_ranking().to_string(index=False))

    print("\nScalePert-Tissue ranking:")
    print(pipeline.tissue_ranking().to_string(index=False))

    print("\nCombined summary table:")
    print(pipeline.summary_table().to_string(index=False))

    pipeline.export("scalepert_results")
    print("\nresults written: scalepert_results_*.csv")


if __name__ == "__main__":
    main()
