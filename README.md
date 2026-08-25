# netComplex

`netComplex` is a standalone Python package for estimating protein-complex
activity from three inputs that researchers can supply themselves:

1. a gene-by-sample expression matrix;
2. a protein-protein interaction (PPI) network; and
3. a table defining protein-complex membership.

The package performs within-sample rank normalisation, propagates gene-level
signal through the PPI using Random Walk with Restart (RWR), and summarises the
propagated values over each complex. It has no bundled biological dataset and
does not require access to any project-specific code or paths.

## When to use netComplex

Use `netComplex` when your scientific question concerns the activity of known
multi-protein complexes rather than individual genes. Typical inputs include
bulk RNA-seq, proteomics, or a processed single-cell expression matrix.

The method is suitable when all three resources use a compatible identifier
namespace, for example HGNC gene symbols in expression, PPI, and complex
membership data.

## Algorithm

Let \(E \in \mathbb{R}^{G \times S}\) be an expression matrix with genes as
rows and samples as columns. For each sample \(s\), netComplex performs:

1. **Within-sample rank normalisation.** Each gene receives the percentile
   rank:

   ```text
   R[g, s] = (rank(E[g, s]) - 0.5) / G
   ```

   This places values in the open interval \((0, 1)\), making samples with
   different expression scales comparable.

2. **Network propagation.** The PPI is converted into an undirected,
   row-normalised adjacency matrix `S`. RWR is iterated until convergence:

   ```text
   X(t + 1) = alpha * R + (1 - alpha) * S * X(t)
   ```

   where \(\alpha\) is the restart probability (`0.3` by default). All genes
   present in the expression matrix are retained as nodes. A gene with no PPI
   edge is an isolated node and keeps its original rank.

3. **Complex scoring.** For complex `c`, the activity score is the mean
   propagated value across its available members:

   ```text
   score[c, s] = mean(X[g, s] for g in members(c) that are present in expression)
   ```

   A complex without any measured member is returned as `NaN`. A partially
   measured complex is scored using its available members; the corresponding
   `coverage` output reports the represented fraction.

By default, retained PPI interactions are unweighted, matching the standard
workflow. Users may instead select a positive numeric confidence column as an
edge-weight parameter; see [Parameters](#parameters).

## Installation

Install the published package:

```bash
pip install netComplex
```

For a reproducible environment, create and activate a virtual environment
before installing:

```bash
python -m venv .venv
# Linux/macOS
source .venv/bin/activate
# Windows PowerShell
# .venv\Scripts\Activate.ps1

python -m pip install --upgrade pip
python -m pip install netComplex
```

`netComplex` supports Python 3.9 or newer and installs `numpy`, `pandas`, and
`scipy` automatically.

## Run the bundled example

After installation, run a complete five-gene, three-sample toy example with no
additional download:

```bash
netcomplex example --output netcomplex-example-results/
```

The command uses a small PPI, an expression matrix, and three example
complexes packaged inside `netComplex`. It writes the same four files produced
by a real analysis:

```text
netcomplex-example-results/
├── complex_scores.csv
├── coverage.csv
├── node_scores.csv
└── run_metadata.json
```

Inspect the score matrix with:

```bash
python -c "import pandas as pd; print(pd.read_csv('netcomplex-example-results/complex_scores.csv', index_col=0))"
```

The bundled `Partially measured example` intentionally contains one gene that
is absent from expression. Its score is calculated from `EGFR`, while its
coverage is `0.5`; this demonstrates how to interpret incomplete complexes.

### Single-cell example

netComplex also scores **each cell independently** when the expression matrix
has cells as columns. The bundled single-cell example contains six cells, five
genes, and realistic zero values from dropout:

```bash
netcomplex example --dataset singlecell --output netcomplex-singlecell-results/
```

`complex_scores.csv` then has complexes as rows and individual cells as
columns. In the bundled data, `cell_01` to `cell_03` have weaker EGFR-pathway
signal than `cell_04` to `cell_06`, giving a small, inspectable per-cell use
case.

For real single-cell data, complete standard quality control, normalisation,
and gene-ID harmonisation before scoring. Use a processed non-negative
gene-by-cell matrix (for example, normalised counts or `log1p`-transformed
counts); do not supply missing values. Because rank normalisation is performed
within each cell, the method is compatible with sparse/dropout-heavy matrices,
but lowly detected complexes should always be interpreted alongside
`coverage.csv`. The package scores cells; cell-type annotation, differential
testing, pseudobulk construction, and batch correction remain downstream
analysis choices.

## Multi-omics extensions

Version `0.3.0` adds two ways to score matched two-omics data, such as RNA and
protein abundance. Both require the two matrices to contain the same genes and
samples. Their row and column order may differ: `netComplex` matches and
reorders identifiers automatically. Missing genes or unmatched samples are an
error rather than being silently discarded.

### 1. RankFusion

`score_rankfusion` rank-normalises each omics layer separately, then uses their
weighted combination as the restart vector on one PPI network:

```text
restart = weight_backbone * rank(backbone)
        + (1 - weight_backbone) * rank(auxiliary)
```

This is the lighter-weight option. The two omics measurements are integrated
before propagation; only one gene-level network is traversed.

```python
from netcomplex import score_rankfusion

rna = pd.read_csv("rna.csv", index_col=0)
protein = pd.read_csv("protein.csv", index_col=0)

result = score_rankfusion(
    backbone_expression=rna,
    auxiliary_expression=protein,
    links=ppi,
    complexes=complexes,
    alpha=0.3,
    weight_backbone=0.5,
)
```

`weight_backbone` must be between `0` and `1`. A value of `0.5` gives the two
layers equal influence; `1.0` uses only the backbone ranks, while retaining the
same PPI calculation.

Command line:

```bash
netcomplex rankfusion \
  --backbone rna.csv \
  --auxiliary protein.csv \
  --ppi ppi.csv \
  --complexes complexes.csv \
  --output rankfusion-results/ \
  --weight-backbone 0.5
```

### 2. Dual-layer PPI

`score_multilayer` creates one PPI layer per omics type. Each gene receives a
same-gene connection between the RNA and protein layers. After RWR, the two
layer scores are collapsed to a gene score before complex averaging:

```text
gene_score = weight_rna * RNA_layer_score
           + (1 - weight_rna) * protein_layer_score
```

`interlayer_weight` controls the strength of the same-gene RNA/protein edge;
it must be greater than zero. The default (`1.0`) gives an inter-layer edge the
same base weight as an intra-layer PPI edge before row normalisation.

```python
from netcomplex import score_multilayer

result = score_multilayer(
    rna_expression=rna,
    protein_expression=protein,
    links=ppi,
    complexes=complexes,
    alpha=0.3,
    weight_rna=0.5,
    interlayer_weight=1.0,
)
```

Command line:

```bash
netcomplex multilayer \
  --rna rna.csv \
  --protein protein.csv \
  --ppi ppi.csv \
  --complexes complexes.csv \
  --output multilayer-results/ \
  --weight-rna 0.5 \
  --interlayer-weight 1.0
```

Both extension APIs return the same `NetComplexResult` fields as `score`:
`complex_scores`, `node_scores`, `coverage`, `iterations`, and `final_delta`.
The command-line metadata records the selected algorithm and its multi-omics
parameters.

Try the bundled protein example with either extension:

```bash
netcomplex example --method rankfusion --output rankfusion-example-results/
netcomplex example --method multilayer --output multilayer-example-results/
```

### Choosing an extension

| Question | Recommended method |
| --- | --- |
| Do you want a simple, transparent weighted integration before one network propagation? | RankFusion (`score_rankfusion`). |
| Do you want RNA and protein signals to exchange information during propagation? | Dual-layer PPI (`score_multilayer`). |
| Do you have unmatched samples or different gene universes? | Neither method yet: match samples and harmonise genes first. |

## Input data requirements

All inputs must use the same gene or protein identifier namespace. The package
does not perform ID conversion, so convert Ensembl IDs, UniProt accessions, or
other identifiers before scoring if necessary.

### Expression matrix

Provide a numeric `pandas.DataFrame` with **genes as the index** and **samples
as columns**. IDs must be unique, values must be finite, and no value may be
missing.

CSV representation:

```text
gene,sample_01,sample_02,sample_03
EGFR,8.2,6.1,7.5
ERBB2,5.4,4.8,6.3
GRB2,3.1,2.9,3.5
```

The first column becomes the DataFrame index when loaded with
`pd.read_csv(..., index_col=0)`.

### PPI network

Provide a DataFrame or CSV with these required columns:

| Column | Meaning |
| --- | --- |
| `protein1` | First interaction endpoint. |
| `protein2` | Second interaction endpoint. |

Additional columns, such as STRING confidence scores, are ignored by default.
Pass their column name as `edge_weight_column` in Python or
`--edge-weight-column` on the command line to use them as positive edge
weights.

```text
protein1,protein2,combined_score
EGFR,ERBB2,0.97
EGFR,GRB2,0.93
```

Self-loops and duplicate edges are removed internally. At least one non-self
PPI edge must overlap the expression gene universe.

### Complex membership table

Provide a DataFrame or CSV containing exactly one row per complex. The `Genes`
column is a semicolon-separated member list.

| Complex | Genes |
| --- | --- |
| EGFR signalling complex | `EGFR;ERBB2;GRB2` |
| Example dimer | `EGFR;ERBB2` |

Complex names must be unique and every complex must list at least one member.
Duplicate gene names within a complex are counted only once.

## Python API

### Basic workflow

```python
import pandas as pd
from netcomplex import score

expression = pd.read_csv("expression.csv", index_col=0)
ppi = pd.read_csv("ppi.csv")
complexes = pd.read_csv("complexes.csv")

result = score(expression, ppi, complexes, alpha=0.3)

# Rows are complexes and columns are samples.
complex_scores = result.complex_scores

# Rows are expression genes and columns are samples.
node_scores = result.node_scores

# Fraction of the listed members available in expression, repeated per sample.
coverage = result.coverage

print(f"RWR iterations: {result.iterations}")
print(f"Final update: {result.final_delta:.3e}")
```

For users who only need the final complex matrix:

```python
from netcomplex import score_matrix

complex_scores = score_matrix(expression, ppi, complexes)
```

Complexes can also be supplied as a Python mapping instead of a DataFrame:

```python
complexes = {
    "EGFR signalling complex": ["EGFR", "ERBB2", "GRB2"],
    "Example dimer": ["EGFR", "ERBB2"],
}
result = score(expression, ppi, complexes)
```

### Parameters

| Parameter | Default | Meaning |
| --- | ---: | --- |
| `alpha` | `0.3` | RWR restart probability. Must be in `(0, 1]`. Larger values preserve more of the original rank; smaller values diffuse further through the PPI. |
| `tolerance` | `1e-4` | Maximum absolute update allowed before RWR is considered converged. Must be positive. |
| `max_iter` | `1000` | Maximum RWR iterations before an error is raised. |
| `rank_method` | `"average"` | Tie rule for within-sample ranks: `"average"`, `"min"`, `"max"`, `"first"`, or `"dense"`. |
| `rank_ascending` | `True` | When `True`, higher expression receives a higher rank. Set `False` to prioritise lower expression. |
| `edge_weight_column` | `None` | Optional PPI column containing finite values greater than zero. When supplied, its values determine relative transition weights. |

RankFusion adds `weight_backbone` (`0.5` by default), the proportion assigned
to the backbone omics rank. Dual-layer scoring adds `weight_rna` (`0.5`) for
collapsing RNA/protein layer scores and `interlayer_weight` (`1.0`) for
same-gene cross-layer edges.

Use the default settings for analyses intended to match the package's standard
rank/RWR/complex-mean workflow. Any non-default setting is written into the
CLI `run_metadata.json`; Python users should record equivalent settings with
the PPI source and version in their analysis metadata.

## Command-line interface

The same calculation is available without writing Python code:

```bash
netcomplex score \
  --expression expression.csv \
  --ppi ppi.csv \
  --complexes complexes.csv \
  --output results/
```

Optional controls:

```bash
netcomplex score \
  --expression expression.csv \
  --ppi ppi.csv \
  --complexes complexes.csv \
  --output results/ \
  --alpha 0.3 \
  --tolerance 1e-4 \
  --max-iter 1000 \
  --rank-method average \
  --edge-weight-column combined_score
```

Use `netcomplex score --help` to view the command contract. The expression CSV
must have gene IDs in its first column; PPI and complex tables must include the
column names described above.

Add `--rank-descending` to make lower expression receive a higher rank. The
same `--rank-method`, `--rank-descending`, and `--edge-weight-column` controls
are available for `netcomplex rankfusion` and `netcomplex multilayer`, together
with their method-specific weight parameters.

### Output files

| File | Rows × columns | Description |
| --- | --- | --- |
| `complex_scores.csv` | complexes × samples | Final NetComplex activity scores. |
| `node_scores.csv` | genes × samples | RWR-propagated gene scores. Useful for diagnostics and downstream interpretation. |
| `coverage.csv` | complexes × samples | Proportion of each complex's listed genes present in expression. |
| `run_metadata.json` | JSON | Algorithm name, parameter values, convergence diagnostics, and output dimensions. |

## Interpretation and reproducibility

- Scores are relative to other genes **within the same sample** because the
  first step is rank normalisation. Interpret cross-sample differences with the
  study design and preprocessing in mind.
- Use a PPI source appropriate to the species and identifier namespace of the
  experiment. Network choice can materially change propagated scores.
- Store the input expression preprocessing, PPI source/version, complex source
  / version, package version, and RWR parameters alongside the output files.
- Review `coverage.csv` before interpreting a complex: a low coverage score
  means its activity was estimated from only a subset of members.
- The package does not infer complexes, perform differential testing, or
  establish causality. Those analyses should be conducted downstream with an
  appropriate statistical design.

## Common errors

| Message | Likely cause | Resolution |
| --- | --- | --- |
| `No expression genes overlap the PPI network` | Different identifier namespaces or species. | Harmonise IDs and verify the PPI source. |
| `No non-self-loop PPI edges remain` | The filtered PPI has only self-loops or no usable interactions. | Supply a PPI with at least one overlap edge between two expression genes. |
| `expression values must be finite` | Missing, infinite, or non-numeric expression values. | Impute or remove missing values before scoring. |
| `Complex names must be unique` | Duplicate rows or duplicate complex labels. | Deduplicate or rename complexes before scoring. |

## Migrating from 0.1.0

PyPI version `0.1.0` documented an experimental implementation using
sample-specific edge reweighting and a coherence multiplier. Version `0.2.0`
standardises the public package on the rank/RWR/complex-mean algorithm above.
This is an intentional algorithm and API change: do not compare scores between
the two versions directly. Version `0.3.0` adds the separate multi-omics APIs
described above. Rerun analyses with the target version and report that version
in the methods section.

## License

Research-only. Contact the author for commercial-use permission.
