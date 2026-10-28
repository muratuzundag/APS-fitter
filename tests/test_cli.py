from pathlib import Path

import numpy as np

from compact_spacing.cli import main
from compact_spacing.io import read_modes

EXAMPLE = Path(__file__).resolve().parents[1] / "examples/example_modes.csv"


def test_csv_groups():
    groups = read_modes(EXAMPLE)
    assert list(groups) == ["l1", "l2"]
    assert len(groups["l1"]["k"]) == 8
    assert "period_error_s" in groups["l1"]


def test_complete_cli(tmp_path):
    assert main([str(EXAMPLE), "--samples", "100", "--output-dir", str(tmp_path),
                 "--format", "png"]) == 0
    assert (tmp_path / "period_spacing_results.csv").is_file()
    for sequence in ["l1", "l2"]:
        assert (tmp_path / f"{sequence}_spacing.png").stat().st_size > 1000
        assert np.load(tmp_path / f"{sequence}_samples.npy").shape == (100,)
