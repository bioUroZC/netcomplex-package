import numpy as np
import pandas as pd
import pytest

from netcomplex import score, score_matrix, score_multilayer, score_rankfusion


@pytest.fixture
def inputs():
    expression = pd.DataFrame(
        {"sample_a": [10.0, 6.0, 2.0, 1.0], "sample_b": [1.0, 5.0, 9.0, 3.0]},
        index=["A", "B", "C", "D"],
    )
    links = pd.DataFrame({"protein1": ["A", "B"], "protein2": ["B", "C"]})
    complexes = pd.DataFrame({"Complex": ["ABC", "CD", "unmeasured"], "Genes": ["A;B;C", "C;D", "X;Y"]})
    return expression, links, complexes


def test_score_returns_expected_shapes_coverage_and_isolated_rank(inputs):
    expression, links, complexes = inputs
    result = score(expression, links, complexes)
    assert result.complex_scores.shape == (3, 2)
    assert result.coverage.loc["ABC", "sample_a"] == 1.0
    assert result.coverage.loc["unmeasured", "sample_a"] == 0.0
    assert result.complex_scores.loc["unmeasured"].isna().all()
    assert result.node_scores.loc["D", "sample_a"] == pytest.approx(0.125)
    assert result.final_delta < 1e-4


def test_score_matrix_matches_full_result(inputs):
    expression, links, complexes = inputs
    pd.testing.assert_frame_equal(score_matrix(expression, links, complexes), score(expression, links, complexes).complex_scores)


def test_expression_must_be_complete_and_numeric(inputs):
    expression, links, complexes = inputs
    expression.loc["A", "sample_a"] = np.nan
    with pytest.raises(ValueError, match="finite"):
        score(expression, links, complexes)


def test_score_accepts_rank_and_ppi_weight_controls(inputs):
    expression, links, complexes = inputs
    links["confidence"] = [0.9, 0.4]
    result = score(
        expression,
        links,
        complexes,
        rank_method="min",
        rank_ascending=False,
        edge_weight_column="confidence",
    )
    assert result.complex_scores.shape == (3, 2)


def test_score_rejects_invalid_ppi_weight_column(inputs):
    expression, links, complexes = inputs
    links["confidence"] = [0.9, 0.0]
    with pytest.raises(ValueError, match="greater than zero"):
        score(expression, links, complexes, edge_weight_column="confidence")


def test_rankfusion_matches_single_omics_for_identical_layers(inputs):
    expression, links, complexes = inputs
    expected = score(expression, links, complexes).complex_scores
    actual = score_rankfusion(expression, expression.copy(), links, complexes).complex_scores
    pd.testing.assert_frame_equal(actual, expected)


def test_multilayer_returns_complex_and_gene_scores(inputs):
    expression, links, complexes = inputs
    result = score_multilayer(expression, expression.copy(), links, complexes)
    assert result.complex_scores.shape == (3, 2)
    assert result.node_scores.shape == (4, 2)
    assert result.coverage.loc["unmeasured", "sample_a"] == 0.0


def test_rankfusion_aligns_equivalent_but_reordered_inputs(inputs):
    expression, links, complexes = inputs
    reordered = expression.loc[["C", "A", "D", "B"], ["sample_b", "sample_a"]]
    result = score_rankfusion(expression, reordered, links, complexes)
    assert list(result.node_scores.index) == ["A", "B", "C", "D"]
    assert list(result.node_scores.columns) == ["sample_a", "sample_b"]
