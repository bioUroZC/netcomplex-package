"""Stable public Python API for netComplex."""

from __future__ import annotations

from dataclasses import dataclass
from collections.abc import Mapping, Sequence

import numpy as np
import pandas as pd

from .core import (
    DEFAULT_ALPHA,
    DEFAULT_MAX_ITER,
    DEFAULT_TOLERANCE,
    aggregate_complexes,
    build_network,
    canonical_links,
    parse_complexes,
    rank_within_samples,
    rwr,
)


@dataclass(frozen=True)
class NetComplexResult:
    """Scores, diagnostics, and convergence information from :func:`score`."""

    complex_scores: pd.DataFrame
    node_scores: pd.DataFrame
    coverage: pd.DataFrame
    iterations: int
    final_delta: float


def _validate_expression(expression: pd.DataFrame) -> pd.DataFrame:
    if not isinstance(expression, pd.DataFrame):
        raise TypeError("expression must be a pandas DataFrame with genes as its index.")
    if expression.empty or expression.shape[1] == 0:
        raise ValueError("expression must contain at least one gene and one sample.")
    if expression.index.has_duplicates or expression.columns.has_duplicates:
        raise ValueError("expression gene and sample identifiers must be unique.")
    if expression.index.isna().any() or expression.columns.isna().any():
        raise ValueError("expression gene and sample identifiers cannot be missing.")
    checked = expression.copy()
    checked.index = checked.index.astype(str).str.strip()
    checked.columns = checked.columns.astype(str).str.strip()
    if (checked.index == "").any() or (checked.columns == "").any():
        raise ValueError("expression gene and sample identifiers cannot be empty.")
    if checked.index.has_duplicates or checked.columns.has_duplicates:
        raise ValueError("Trimming identifiers produced duplicates; make identifiers unique first.")
    try:
        checked = checked.astype(float)
    except (TypeError, ValueError) as exc:
        raise ValueError("expression values must be numeric.") from exc
    if not np.isfinite(checked.to_numpy()).all():
        raise ValueError("expression values must be finite and contain no missing values.")
    return checked


def score(
    expression: pd.DataFrame,
    links: pd.DataFrame,
    complexes: pd.DataFrame | Mapping[str, Sequence[str]],
    *,
    alpha: float = DEFAULT_ALPHA,
    tolerance: float = DEFAULT_TOLERANCE,
    max_iter: int = DEFAULT_MAX_ITER,
    rank_method: str = "average",
    rank_ascending: bool = True,
    edge_weight_column: str | None = None,
) -> NetComplexResult:
    """Score all complexes across all samples in user-supplied data.

    ``expression`` is a numeric genes-by-samples DataFrame. ``links`` has
    ``protein1`` and ``protein2`` columns. ``complexes`` is either a DataFrame
    containing ``Complex`` and semicolon-delimited ``Genes`` columns or a
    mapping from complex name to its ordered member genes. ``rank_method`` and
    ``rank_ascending`` control within-sample ranking. Set
    ``edge_weight_column`` to a positive numeric PPI column name to use a
    weighted interaction network.
    """
    checked_expression = _validate_expression(expression)
    checked_links = canonical_links(links, edge_weight_column=edge_weight_column)
    parsed_complexes = parse_complexes(complexes)
    nodes, smoothing, degree = build_network(checked_links, checked_expression.index)
    ranks = rank_within_samples(checked_expression, nodes, method=rank_method, ascending=rank_ascending)
    node_scores, iterations, final_delta = rwr(
        ranks, smoothing, degree, alpha=alpha, tolerance=tolerance, max_iter=max_iter
    )
    complex_scores, coverage = aggregate_complexes(node_scores, parsed_complexes)
    return NetComplexResult(complex_scores, node_scores, coverage, iterations, final_delta)


def score_matrix(
    expression: pd.DataFrame,
    links: pd.DataFrame,
    complexes: pd.DataFrame | Mapping[str, Sequence[str]],
    *,
    alpha: float = DEFAULT_ALPHA,
    tolerance: float = DEFAULT_TOLERANCE,
    max_iter: int = DEFAULT_MAX_ITER,
    rank_method: str = "average",
    rank_ascending: bool = True,
    edge_weight_column: str | None = None,
) -> pd.DataFrame:
    """Return only the complexes-by-samples activity matrix."""
    return score(
        expression, links, complexes, alpha=alpha, tolerance=tolerance, max_iter=max_iter,
        rank_method=rank_method, rank_ascending=rank_ascending, edge_weight_column=edge_weight_column,
    ).complex_scores
