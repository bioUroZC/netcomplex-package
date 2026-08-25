"""Public multi-omics extensions for netComplex."""

from __future__ import annotations

from collections.abc import Mapping, Sequence

import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix

from .api import NetComplexResult, _validate_expression
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


def _align_omics(
    primary: pd.DataFrame, secondary: pd.DataFrame, *, primary_name: str, secondary_name: str
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Validate two matrices and align secondary IDs to the primary matrix."""
    checked_primary = _validate_expression(primary)
    checked_secondary = _validate_expression(secondary)
    primary_genes, secondary_genes = set(checked_primary.index), set(checked_secondary.index)
    primary_samples, secondary_samples = set(checked_primary.columns), set(checked_secondary.columns)
    if primary_genes != secondary_genes:
        raise ValueError(
            f"{primary_name} and {secondary_name} must contain the same gene identifiers; "
            f"only in {primary_name}: {sorted(primary_genes - secondary_genes)[:5]}, "
            f"only in {secondary_name}: {sorted(secondary_genes - primary_genes)[:5]}."
        )
    if primary_samples != secondary_samples:
        raise ValueError(
            f"{primary_name} and {secondary_name} must contain the same sample identifiers; "
            f"only in {primary_name}: {sorted(primary_samples - secondary_samples)[:5]}, "
            f"only in {secondary_name}: {sorted(secondary_samples - primary_samples)[:5]}."
        )
    return checked_primary, checked_secondary.loc[checked_primary.index, checked_primary.columns]


def _check_unit_interval(value: float, name: str) -> None:
    if not 0.0 <= value <= 1.0:
        raise ValueError(f"{name} must lie in [0, 1].")


def score_rankfusion(
    backbone_expression: pd.DataFrame,
    auxiliary_expression: pd.DataFrame,
    links: pd.DataFrame,
    complexes: pd.DataFrame | Mapping[str, Sequence[str]],
    *,
    alpha: float = DEFAULT_ALPHA,
    weight_backbone: float = 0.5,
    tolerance: float = DEFAULT_TOLERANCE,
    max_iter: int = DEFAULT_MAX_ITER,
    rank_method: str = "average",
    rank_ascending: bool = True,
    edge_weight_column: str | None = None,
) -> NetComplexResult:
    """Score complexes using RankFusion of two matched omics layers.

    The two matrices must represent the same genes and samples; their row and
    column order may differ and is aligned automatically. Each layer is ranked
    within sample. The RWR restart matrix is then
    ``weight_backbone * ranks_backbone + (1 - weight_backbone) * ranks_auxiliary``
    on one fixed, unweighted PPI network.
    """
    _check_unit_interval(weight_backbone, "weight_backbone")
    backbone, auxiliary = _align_omics(
        backbone_expression,
        auxiliary_expression,
        primary_name="backbone_expression",
        secondary_name="auxiliary_expression",
    )
    checked_links = canonical_links(links, edge_weight_column=edge_weight_column)
    parsed_complexes = parse_complexes(complexes)
    nodes, smoothing, degree = build_network(checked_links, backbone.index)
    backbone_ranks = rank_within_samples(backbone, nodes, method=rank_method, ascending=rank_ascending)
    auxiliary_ranks = rank_within_samples(auxiliary, nodes, method=rank_method, ascending=rank_ascending)
    restart = weight_backbone * backbone_ranks + (1.0 - weight_backbone) * auxiliary_ranks
    node_scores, iterations, final_delta = rwr(
        restart, smoothing, degree, alpha=alpha, tolerance=tolerance, max_iter=max_iter
    )
    complex_scores, coverage = aggregate_complexes(node_scores, parsed_complexes)
    return NetComplexResult(complex_scores, node_scores, coverage, iterations, final_delta)


def _build_multilayer_network(
    links: pd.DataFrame, genes: pd.Index, interlayer_weight: float
) -> tuple[pd.Index, csr_matrix, np.ndarray]:
    """Build a row-normalised RNA/protein PPI graph with inter-layer edges."""
    if interlayer_weight <= 0:
        raise ValueError("interlayer_weight must be positive.")
    nodes, _, _ = build_network(links, genes)
    node_position = {gene: position for position, gene in enumerate(nodes)}
    edge_table = links[links.protein1.isin(nodes) & links.protein2.isin(nodes)].copy()
    edge_table = edge_table[edge_table.protein1 != edge_table.protein2].drop_duplicates()
    n_nodes = len(nodes)
    row: list[int] = []
    col: list[int] = []
    data: list[float] = []

    def add_undirected(left: int, right: int, weight: float) -> None:
        row.extend((left, right))
        col.extend((right, left))
        data.extend((weight, weight))

    edge_weights = edge_table["_edge_weight"].to_numpy(dtype=float) if "_edge_weight" in edge_table else np.ones(len(edge_table))
    for left_gene, right_gene, edge_weight in zip(edge_table.protein1.to_numpy(), edge_table.protein2.to_numpy(), edge_weights):
        left, right = node_position[left_gene], node_position[right_gene]
        add_undirected(left, right, float(edge_weight))
        add_undirected(left + n_nodes, right + n_nodes, float(edge_weight))
    for position in range(n_nodes):
        add_undirected(position, position + n_nodes, interlayer_weight)

    adjacency = csr_matrix((data, (row, col)), shape=(2 * n_nodes, 2 * n_nodes))
    degree = np.asarray(adjacency.sum(axis=1)).ravel()
    inverse_degree = np.divide(1.0, degree, out=np.zeros_like(degree), where=degree > 0)
    return nodes, adjacency.multiply(inverse_degree[:, None]).tocsr(), degree


def score_multilayer(
    rna_expression: pd.DataFrame,
    protein_expression: pd.DataFrame,
    links: pd.DataFrame,
    complexes: pd.DataFrame | Mapping[str, Sequence[str]],
    *,
    alpha: float = DEFAULT_ALPHA,
    weight_rna: float = 0.5,
    interlayer_weight: float = 1.0,
    tolerance: float = DEFAULT_TOLERANCE,
    max_iter: int = DEFAULT_MAX_ITER,
    rank_method: str = "average",
    rank_ascending: bool = True,
    edge_weight_column: str | None = None,
) -> NetComplexResult:
    """Score complexes with an RNA/protein dual-layer PPI network.

    Matched RNA and protein matrices each occupy one copy of the PPI. Every
    gene also receives an RNA-to-protein edge with ``interlayer_weight``. After
    RWR, the two layer scores are collapsed as
    ``weight_rna * RNA + (1 - weight_rna) * protein`` before complex scoring.
    """
    _check_unit_interval(weight_rna, "weight_rna")
    rna, protein = _align_omics(
        rna_expression, protein_expression, primary_name="rna_expression", secondary_name="protein_expression"
    )
    checked_links = canonical_links(links, edge_weight_column=edge_weight_column)
    parsed_complexes = parse_complexes(complexes)
    nodes, smoothing, degree = _build_multilayer_network(checked_links, rna.index, interlayer_weight)
    rna_ranks = rank_within_samples(rna, nodes, method=rank_method, ascending=rank_ascending)
    protein_ranks = rank_within_samples(protein, nodes, method=rank_method, ascending=rank_ascending)
    restart = pd.concat((rna_ranks, protein_ranks), axis=0)
    restart.index = [f"RNA::{gene}" for gene in nodes] + [f"PROTEIN::{gene}" for gene in nodes]
    multilayer_scores, iterations, final_delta = rwr(
        restart, smoothing, degree, alpha=alpha, tolerance=tolerance, max_iter=max_iter
    )
    n_nodes = len(nodes)
    collapsed = (
        weight_rna * multilayer_scores.iloc[:n_nodes].to_numpy(dtype=float)
        + (1.0 - weight_rna) * multilayer_scores.iloc[n_nodes:].to_numpy(dtype=float)
    )
    node_scores = pd.DataFrame(collapsed, index=nodes, columns=rna.columns)
    complex_scores, coverage = aggregate_complexes(node_scores, parsed_complexes)
    return NetComplexResult(complex_scores, node_scores, coverage, iterations, final_delta)
