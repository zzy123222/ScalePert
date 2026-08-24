# ScalePert development notes

## Layout

- `src/scalepert/` — package source (Cell and Tissue modules, pipeline)
- `tests/` — pytest suite covering displacement math, sign conventions, propagation
- `examples/` — runnable end-to-end scripts
- `scripts/build_example.py` — regenerates `data/scalepert_example.h5ad`
- `data/` — bundled example atlas

## Local workflow

```bash
python -m venv .venv && source .venv/bin/activate
pip install -e ".[dev]"
pytest -q
python examples/run_pipeline.py
```

To rebuild the example dataset:

```bash
python scripts/build_example.py data/scalepert_example.h5ad
```

## Conventions

- No inline comments; names and docstrings carry the semantics.
- Default parameters match the manuscript: k=30, s=1.0, beta=0.3, 50 PCs.
- Scores: signed Wasserstein-1 shifts (negative = program suppression);
  tissue layer uses non-negative suppression magnitude.
