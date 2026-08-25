import json

import pandas as pd

from netcomplex.cli import main


def test_cli_writes_result_bundle(tmp_path):
    pd.DataFrame({"sample": [5.0, 2.0]}, index=["A", "B"]).to_csv(tmp_path / "expression.csv")
    pd.DataFrame({"protein1": ["A"], "protein2": ["B"], "confidence": [0.8]}).to_csv(tmp_path / "ppi.csv", index=False)
    pd.DataFrame({"Complex": ["AB"], "Genes": ["A;B"]}).to_csv(tmp_path / "complexes.csv", index=False)
    output = tmp_path / "results"
    assert main(["score", "--expression", str(tmp_path / "expression.csv"), "--ppi", str(tmp_path / "ppi.csv"), "--complexes", str(tmp_path / "complexes.csv"), "--output", str(output)]) == 0
    metadata = json.loads((output / "run_metadata.json").read_text(encoding="utf-8"))
    assert metadata["algorithm"] == "rank-rwr-complex-mean"
    assert (output / "complex_scores.csv").is_file()


def test_bundled_example_writes_result_bundle(tmp_path):
    output = tmp_path / "example-results"
    assert main(["example", "--output", str(output)]) == 0
    scores = pd.read_csv(output / "complex_scores.csv", index_col=0)
    assert scores.shape == (3, 3)
    assert (output / "coverage.csv").is_file()


def test_bundled_multilayer_example_writes_result_bundle(tmp_path):
    output = tmp_path / "multilayer-example-results"
    assert main(["example", "--method", "multilayer", "--output", str(output)]) == 0
    metadata = json.loads((output / "run_metadata.json").read_text(encoding="utf-8"))
    assert metadata["algorithm"] == "dual-layer-ppi-multiomics"


def test_bundled_singlecell_example_writes_per_cell_scores(tmp_path):
    output = tmp_path / "singlecell-example-results"
    assert main(["example", "--dataset", "singlecell", "--output", str(output)]) == 0
    scores = pd.read_csv(output / "complex_scores.csv", index_col=0)
    assert scores.shape == (3, 6)
    assert list(scores.columns) == ["cell_01", "cell_02", "cell_03", "cell_04", "cell_05", "cell_06"]


def test_rankfusion_command_writes_metadata(tmp_path):
    pd.DataFrame({"sample": [5.0, 2.0]}, index=["A", "B"]).to_csv(tmp_path / "rna.csv")
    pd.DataFrame({"sample": [3.0, 8.0]}, index=["A", "B"]).to_csv(tmp_path / "protein.csv")
    pd.DataFrame({"protein1": ["A"], "protein2": ["B"], "confidence": [0.8]}).to_csv(tmp_path / "ppi.csv", index=False)
    pd.DataFrame({"Complex": ["AB"], "Genes": ["A;B"]}).to_csv(tmp_path / "complexes.csv", index=False)
    output = tmp_path / "rankfusion-results"
    assert main(["rankfusion", "--backbone", str(tmp_path / "rna.csv"), "--auxiliary", str(tmp_path / "protein.csv"), "--ppi", str(tmp_path / "ppi.csv"), "--complexes", str(tmp_path / "complexes.csv"), "--output", str(output), "--rank-method", "min", "--rank-descending", "--edge-weight-column", "confidence"]) == 0
    metadata = json.loads((output / "run_metadata.json").read_text(encoding="utf-8"))
    assert metadata["algorithm"] == "rankfusion"
    assert metadata["weight_backbone"] == 0.5
    assert metadata["rank_method"] == "min"
    assert metadata["rank_ascending"] is False
    assert metadata["edge_weight_column"] == "confidence"
