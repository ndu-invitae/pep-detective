import pandas as pd
import pytest

from pep_detective.src.ancova import AncovaResult, CovarProcessor
from tests.pep_detective.conftest import TEST_DATA_PATH


@pytest.mark.parametrize(
    "sample, tsv_path, pval, expected",
    [
        [
            "test_sample1",
            TEST_DATA_PATH / "test_sample1.tsv",
            0.05,
            AncovaResult(
                sample_id="test_sample1",
                ph_covar=True,
                p_enhancer=0.035128305312971356,
                p_suppressor=0.9648716946870287,
            ),
        ],
        [
            "test_sample2",
            TEST_DATA_PATH / "test_sample2.tsv",
            0.05,
            AncovaResult(
                sample_id="test_sample2",
                ph_covar=False,
                p_enhancer=0.04831824569789305,
                p_suppressor=0.9516817543021069,
            ),
        ],
    ],
)
def test_CovarProcessor(sample, tsv_path, pval, expected):
    """Verify proper ANCOVA and t test results."""
    analysis_results = CovarProcessor(sample, tsv_path, pval).t_test()
    assert analysis_results == expected
