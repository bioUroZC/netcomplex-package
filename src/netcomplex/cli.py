"""Command-line access to the standalone netComplex package."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from importlib import resources

import pandas as pd

from .api import score
from .multiomics import score_multilayer, score_rankfusion


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="netcomplex", description="Score protein complexes with rank/RWR propagation.")
    commands = parser.add_subparsers(dest="command", required=True)
    command = commands.add_parser("score", help="score CSV expression, PPI, and complex-membership inputs")
    command.add_argument("--expression", type=Path, required=True, help="genes-by-samples CSV; first column is the gene ID")
    command.add_argument("--ppi", type=Path, required=True, help="PPI CSV with protein1 and protein2 columns")
    command.add_argument("--complexes", type=Path, required=True, help="complex CSV with Complex and Genes columns")
    command.add_argument("--output", type=Path, required=True, help="directory to receive result files")
    command.add_argument("--alpha", type=float, default=0.3, help="restart probability (default: 0.3)")
    command.add_argument("--tolerance", type=float, default=1e-4, help="convergence tolerance (default: 1e-4)")
    command.add_argument("--max-iter", type=int, default=1000, help="maximum RWR iterations (default: 1000)")
    command.add_argument("--rank-method", choices=("average", "min", "max", "first", "dense"), default="average", help="within-sample tie-ranking method (default: average)")
    command.add_argument("--rank-descending", action="store_false", dest="rank_ascending", help="assign higher ranks to lower expression values")
    command.set_defaults(rank_ascending=True)
    command.add_argument("--edge-weight-column", default=None, help="optional positive numeric PPI column to use as edge weights")
    fusion = commands.add_parser("rankfusion", help="score matched omics by rank fusion before RWR")
    fusion.add_argument("--backbone", type=Path, required=True, help="backbone omics genes-by-samples CSV")
    fusion.add_argument("--auxiliary", type=Path, required=True, help="auxiliary omics genes-by-samples CSV")
    fusion.add_argument("--ppi", type=Path, required=True, help="PPI CSV with protein1 and protein2 columns")
    fusion.add_argument("--complexes", type=Path, required=True, help="complex CSV with Complex and Genes columns")
    fusion.add_argument("--output", type=Path, required=True, help="directory to receive result files")
    fusion.add_argument("--alpha", type=float, default=0.3, help="restart probability (default: 0.3)")
    fusion.add_argument("--weight-backbone", type=float, default=0.5, help="backbone rank weight (default: 0.5)")
    fusion.add_argument("--tolerance", type=float, default=1e-4, help="convergence tolerance (default: 1e-4)")
    fusion.add_argument("--max-iter", type=int, default=1000, help="maximum RWR iterations (default: 1000)")
    fusion.add_argument("--rank-method", choices=("average", "min", "max", "first", "dense"), default="average", help="within-sample tie-ranking method (default: average)")
    fusion.add_argument("--rank-descending", action="store_false", dest="rank_ascending", help="assign higher ranks to lower expression values")
    fusion.set_defaults(rank_ascending=True)
    fusion.add_argument("--edge-weight-column", default=None, help="optional positive numeric PPI column to use as edge weights")
    multilayer = commands.add_parser("multilayer", help="score matched RNA/protein data on a dual-layer PPI")
    multilayer.add_argument("--rna", type=Path, required=True, help="RNA genes-by-samples CSV")
    multilayer.add_argument("--protein", type=Path, required=True, help="protein genes-by-samples CSV")
    multilayer.add_argument("--ppi", type=Path, required=True, help="PPI CSV with protein1 and protein2 columns")
    multilayer.add_argument("--complexes", type=Path, required=True, help="complex CSV with Complex and Genes columns")
    multilayer.add_argument("--output", type=Path, required=True, help="directory to receive result files")
    multilayer.add_argument("--alpha", type=float, default=0.3, help="restart probability (default: 0.3)")
    multilayer.add_argument("--weight-rna", type=float, default=0.5, help="RNA collapse weight (default: 0.5)")
    multilayer.add_argument("--interlayer-weight", type=float, default=1.0, help="same-gene cross-layer edge weight (default: 1.0)")
    multilayer.add_argument("--tolerance", type=float, default=1e-4, help="convergence tolerance (default: 1e-4)")
    multilayer.add_argument("--max-iter", type=int, default=1000, help="maximum RWR iterations (default: 1000)")
    multilayer.add_argument("--rank-method", choices=("average", "min", "max", "first", "dense"), default="average", help="within-sample tie-ranking method (default: average)")
    multilayer.add_argument("--rank-descending", action="store_false", dest="rank_ascending", help="assign higher ranks to lower expression values")
    multilayer.set_defaults(rank_ascending=True)
    multilayer.add_argument("--edge-weight-column", default=None, help="optional positive numeric PPI column to use as edge weights")
    example = commands.add_parser("example", help="run the bundled toy dataset")
    example.add_argument("--output", type=Path, required=True, help="directory to receive example result files")
    example.add_argument("--dataset", choices=("bulk", "singlecell"), default="bulk", help="bundled dataset to use (default: bulk)")
    example.add_argument("--method", choices=("single", "rankfusion", "multilayer"), default="single", help="scoring method to demonstrate (default: single)")
    return parser


def _write_results(result, output: Path, *, algorithm: str, alpha: float, tolerance: float, max_iter: int, **parameters: object) -> None:
    """Write a standard result bundle for either user or bundled inputs."""
    output.mkdir(parents=True, exist_ok=True)
    result.complex_scores.to_csv(output / "complex_scores.csv")
    result.node_scores.to_csv(output / "node_scores.csv")
    result.coverage.to_csv(output / "coverage.csv")
    (output / "run_metadata.json").write_text(
        json.dumps({
            "algorithm": algorithm,
            "alpha": alpha,
            "tolerance": tolerance,
            "max_iter": max_iter,
            "iterations": result.iterations,
            "final_delta": result.final_delta,
            "n_genes": int(result.node_scores.shape[0]),
            "n_samples": int(result.node_scores.shape[1]),
            "n_complexes": int(result.complex_scores.shape[0]),
            **parameters,
        }, indent=2) + "\n",
        encoding="utf-8",
    )
    print(f"Scored {result.complex_scores.shape[0]} complexes across {result.node_scores.shape[1]} samples.")
    print(f"RWR converged in {result.iterations} iterations (delta={result.final_delta:.3e}).")


def main(argv: list[str] | None = None) -> int:
    """Execute the CLI and return its status code."""
    args = _parser().parse_args(argv)
    if args.command == "score":
        expression = pd.read_csv(args.expression, index_col=0)
        links = pd.read_csv(args.ppi)
        complexes = pd.read_csv(args.complexes)
        result = score(expression, links, complexes, alpha=args.alpha, tolerance=args.tolerance, max_iter=args.max_iter, rank_method=args.rank_method, rank_ascending=args.rank_ascending, edge_weight_column=args.edge_weight_column)
        _write_results(result, args.output, algorithm="rank-rwr-complex-mean", alpha=args.alpha, tolerance=args.tolerance, max_iter=args.max_iter, rank_method=args.rank_method, rank_ascending=args.rank_ascending, edge_weight_column=args.edge_weight_column)
        return 0
    if args.command == "rankfusion":
        backbone = pd.read_csv(args.backbone, index_col=0)
        auxiliary = pd.read_csv(args.auxiliary, index_col=0)
        links = pd.read_csv(args.ppi)
        complexes = pd.read_csv(args.complexes)
        result = score_rankfusion(backbone, auxiliary, links, complexes, alpha=args.alpha, weight_backbone=args.weight_backbone, tolerance=args.tolerance, max_iter=args.max_iter, rank_method=args.rank_method, rank_ascending=args.rank_ascending, edge_weight_column=args.edge_weight_column)
        _write_results(result, args.output, algorithm="rankfusion", alpha=args.alpha, tolerance=args.tolerance, max_iter=args.max_iter, weight_backbone=args.weight_backbone, rank_method=args.rank_method, rank_ascending=args.rank_ascending, edge_weight_column=args.edge_weight_column)
        return 0
    if args.command == "multilayer":
        rna = pd.read_csv(args.rna, index_col=0)
        protein = pd.read_csv(args.protein, index_col=0)
        links = pd.read_csv(args.ppi)
        complexes = pd.read_csv(args.complexes)
        result = score_multilayer(rna, protein, links, complexes, alpha=args.alpha, weight_rna=args.weight_rna, interlayer_weight=args.interlayer_weight, tolerance=args.tolerance, max_iter=args.max_iter, rank_method=args.rank_method, rank_ascending=args.rank_ascending, edge_weight_column=args.edge_weight_column)
        _write_results(result, args.output, algorithm="dual-layer-ppi-multiomics", alpha=args.alpha, tolerance=args.tolerance, max_iter=args.max_iter, weight_rna=args.weight_rna, interlayer_weight=args.interlayer_weight, rank_method=args.rank_method, rank_ascending=args.rank_ascending, edge_weight_column=args.edge_weight_column)
        return 0

    with resources.as_file(resources.files("netcomplex").joinpath("examples")) as example_dir:
        expression_file = "singlecell_expression.csv" if args.dataset == "singlecell" else "expression.csv"
        expression = pd.read_csv(example_dir / expression_file, index_col=0)
        links = pd.read_csv(example_dir / "ppi.csv")
        complexes = pd.read_csv(example_dir / "complexes.csv")
        protein = pd.read_csv(example_dir / "protein_expression.csv", index_col=0) if args.dataset == "bulk" else None
    if args.dataset == "singlecell" and args.method != "single":
        raise ValueError("The bundled single-cell dataset demonstrates single-omics scoring; use --method single.")
    if args.method == "rankfusion":
        result = score_rankfusion(expression, protein, links, complexes)
        _write_results(result, args.output, algorithm="rankfusion", alpha=0.3, tolerance=1e-4, max_iter=1000, weight_backbone=0.5)
    elif args.method == "multilayer":
        result = score_multilayer(expression, protein, links, complexes)
        _write_results(result, args.output, algorithm="dual-layer-ppi-multiomics", alpha=0.3, tolerance=1e-4, max_iter=1000, weight_rna=0.5, interlayer_weight=1.0)
    else:
        result = score(expression, links, complexes)
        _write_results(result, args.output, algorithm="rank-rwr-complex-mean", alpha=0.3, tolerance=1e-4, max_iter=1000)
    print(f"Used the bundled {args.dataset} example dataset ({args.method}).")
    return 0
