"""Numerical implementation of the rank/RWR/complex-mean algorithm."""

from __future__ import annotations

from collections.abc import Mapping, Sequence

import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix

DEFAULT_ALPHA = 0.3
DEFAULT_TOLERANCE = 1e-4
DEFAULT_MAX_ITER = 1000
RANK_METHODS = ("average", "min", "max", "first", "dense")


def canonical_links(links: pd.DataFrame, *, edge_weight_column: str | None = None) -> pd.DataFrame:
    """Validate PPI input and retain its endpoint columns.

    Edges are unweighted by default. Set ``edge_weight_column`` to use a
    positive numeric PPI column as an edge weight during row normalisation.
    """
    if not isinstance(links, pd.DataFrame):
        raise TypeError("links must be a pandas DataFrame.")
    required = {"protein1", "protein2"}
    if edge_weight_column is not None:
        required.add(edge_weight_column)
    missing = required.difference(links.columns)
    if missing:
        raise ValueError(f"PPI links missing required column(s): {sorted(missing)}")
    columns = ["protein1", "protein2"] + ([edge_weight_column] if edge_weight_column else [])
    result = links.loc[:, columns].copy()
    if result[["protein1", "protein2"]].isna().any().any():
        raise ValueError("PPI links contain missing protein identifiers.")
    for column in ("protein1", "protein2"):
        result[column] = result[column].astype(str).str.strip()
    if (result[["protein1", "protein2"]] == "").any().any():
        raise ValueError("PPI links contain empty protein identifiers.")
    if edge_weight_column is not None:
        if result[edge_weight_column].isna().any():
            raise ValueError(f"PPI edge weight column {edge_weight_column!r} contains missing values.")
        try:
            weights = pd.to_numeric(result[edge_weight_column], errors="raise").to_numpy(dtype=float)
        except (TypeError, ValueError) as exc:
            raise ValueError(f"PPI edge weight column {edge_weight_column!r} must be numeric.") from exc
        if not np.isfinite(weights).all() or (weights <= 0).any():
            raise ValueError(f"PPI edge weight column {edge_weight_column!r} must contain finite values greater than zero.")
        result = result.rename(columns={edge_weight_column: "_edge_weight"})
    return result


def build_network(links: pd.DataFrame, genes: pd.Index) -> tuple[pd.Index, csr_matrix, np.ndarray]:
    """Construct a row-normalised graph across all supplied expression genes."""
    nodes = pd.Index(sorted(genes))
    node_position = {gene: position for position, gene in enumerate(nodes)}
    linked_genes = pd.Index(pd.unique(links[["protein1", "protein2"]].to_numpy().ravel()))
    overlap = nodes.intersection(linked_genes)
    if overlap.empty:
        raise ValueError("No expression genes overlap the PPI network.")

    edges = links[links.protein1.isin(nodes) & links.protein2.isin(nodes)].copy()
    edges = edges[edges.protein1 != edges.protein2].drop_duplicates()
    if edges.empty:
        raise ValueError("No non-self-loop PPI edges remain after expression filtering.")
    row = np.concatenate((edges.protein1.map(node_position).to_numpy(), edges.protein2.map(node_position).to_numpy()))
    col = np.concatenate((edges.protein2.map(node_position).to_numpy(), edges.protein1.map(node_position).to_numpy()))
    edge_weights = edges["_edge_weight"].to_numpy(dtype=float) if "_edge_weight" in edges else np.ones(len(edges))
    adjacency = csr_matrix((np.concatenate((edge_weights, edge_weights)), (row, col)), shape=(len(nodes), len(nodes)))
    degree = np.asarray(adjacency.sum(axis=1)).ravel()
    inverse_degree = np.divide(1.0, degree, out=np.zeros_like(degree), where=degree > 0)
    return nodes, adjacency.multiply(inverse_degree[:, None]).tocsr(), degree


def rank_within_samples(
    expression: pd.DataFrame,
    nodes: pd.Index,
    *,
    method: str = "average",
    ascending: bool = True,
) -> pd.DataFrame:
    """Return percentile gene ranks in the open interval ``(0, 1)``."""
    if method not in RANK_METHODS:
        raise ValueError(f"rank_method must be one of {list(RANK_METHODS)}, got {method!r}.")
    values = expression.loc[nodes]
    ranks = values.rank(axis=0, method=method, ascending=ascending)
    return (ranks - 0.5) / len(nodes)


def rwr(
    ranks: pd.DataFrame,
    smoothing: csr_matrix,
    degree: np.ndarray,
    *,
    alpha: float,
    tolerance: float,
    max_iter: int,
) -> tuple[pd.DataFrame, int, float]:
    """Diffuse ranks with restart, preserving ranks of isolated genes."""
    if not 0.0 < alpha <= 1.0:
        raise ValueError("alpha must lie in (0, 1].")
    if tolerance <= 0:
        raise ValueError("tolerance must be positive.")
    if max_iter < 1:
        raise ValueError("max_iter must be at least one.")
    original = ranks.to_numpy(dtype=float)
    propagated = original.copy()
    isolated = degree == 0
    for iteration in range(1, max_iter + 1):
        updated = alpha * original + (1.0 - alpha) * (smoothing @ propagated)
        updated[isolated, :] = original[isolated, :]
        delta = float(np.max(np.abs(updated - propagated)))
        propagated = updated
        if delta < tolerance:
            return pd.DataFrame(propagated, index=ranks.index, columns=ranks.columns), iteration, delta
    raise RuntimeError(f"RWR did not converge within {max_iter} iterations.")


def parse_complexes(complexes: pd.DataFrame | Mapping[str, Sequence[str]]) -> dict[str, list[str]]:
    """Convert standard tabular or mapping membership data to member lists."""
    if isinstance(complexes, pd.DataFrame):
        required = {"Complex", "Genes"}
        missing = required.difference(complexes.columns)
        if missing:
            raise ValueError(f"Complex table missing required column(s): {sorted(missing)}")
        records = ((str(row.Complex), str(row.Genes).split(";")) for row in complexes.itertuples(index=False))
    elif isinstance(complexes, Mapping):
        records = ((str(name), members) for name, members in complexes.items())
    else:
        raise TypeError("complexes must be a DataFrame or mapping of complex names to gene sequences.")

    parsed: dict[str, list[str]] = {}
    for name, members in records:
        if not name or name.lower() == "nan":
            raise ValueError("Complex names must be non-empty strings.")
        if name in parsed:
            raise ValueError(f"Complex names must be unique; found duplicate {name!r}.")
        unique_members = list(dict.fromkeys(str(gene).strip() for gene in members if str(gene).strip()))
        if not unique_members:
            raise ValueError(f"Complex {name!r} has no gene members.")
        parsed[name] = unique_members
    if not parsed:
        raise ValueError("At least one complex is required.")
    return parsed


def aggregate_complexes(node_scores: pd.DataFrame, complexes: Mapping[str, Sequence[str]]) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Compute complex means and repeatable per-sample member coverage."""
    available_genes = set(node_scores.index)
    scores: dict[str, pd.Series] = {}
    coverage: dict[str, float] = {}
    for name, members in complexes.items():
        available = [gene for gene in members if gene in available_genes]
        coverage[name] = len(available) / len(members)
        scores[name] = node_scores.loc[available].mean(axis=0) if available else pd.Series(np.nan, index=node_scores.columns)
    score_matrix = pd.DataFrame(scores).T
    coverage_matrix = pd.DataFrame(
        np.repeat(np.asarray(list(coverage.values()))[:, None], node_scores.shape[1], axis=1),
        index=list(coverage),
        columns=node_scores.columns,
    )
    score_matrix.index.name = coverage_matrix.index.name = "Complex"
    return score_matrix, coverage_matrix
