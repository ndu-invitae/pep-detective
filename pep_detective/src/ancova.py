from __future__ import annotations

import logging
from dataclasses import dataclass

import numpy as np
import pandas as pd
from pingouin import ancova
from statsmodels.formula.api import ols
from statsmodels.stats.weightstats import ttest_ind

log = logging.getLogger(__name__)


@dataclass
class AncovaResult:
    """Dataclass for storing ANCOVA analysis results."""

    sample_id: str
    ph_covar: bool
    p_enhancer: float
    p_suppressor: float


class CovarProcessor:
    """
    Main class for processing activity data.

    input_path: Input tsv file that contains one dependent variable (dv) column, one between-subjects factor column with two levels (bf), and one covariate (covar) column
    max_covar_p_value: maximal allowed P value for rejecting no significant correlation hypothesis (default 0.05)
    """

    def __init__(
        self,
        sample_id: str,
        input_path: str,
        max_covar_p_value: float = 0.05,
    ):

        self.sample_id = sample_id
        self.input_path = input_path
        self.max_covar_p_value = max_covar_p_value

        self.open_tsv()

    def open_tsv(self):
        """load in tsv file"""
        try:
            self.df_raw = pd.read_table(self.input_path)
        except Exception as e:
            log.error(
                'Failed to load TSV file: "%s". Error: %s',
                self.sample_id,
                str(e),
            )
            raise

    def ancova_analysis(self) -> bool:
        """Run ANCOVA and determine if pH is a significant covariate"""

        covar_stats = ancova(data=self.df_raw, dv="activity", between="treatment", covar="ph")
        p_ph = covar_stats.query('Source == "ph"')["p-unc"].values[0]

        if np.isnan(p_ph):
            log.warning(
                "Failed to get p value for pH correlation",
                self.sample_id,
            )
        log.info("ANCOVA results", self.sample_id, covar_stats.to_markdown())
        return p_ph < self.max_covar_p_value

    def t_test(self) -> AncovaResult:
        """Execute t test on activities with/ without peptide treatment and return"""

        if (
            significant_ph := self.ancova_analysis()
        ):  # perform linear fitting if pH is a significant covariate otherwise use mean for residue calculaition
            lm = ols("activity ~ ph", self.df_raw).fit()
            y_hat = lm.predict(self.df_raw["ph"])
            self.df_raw["residual"] = self.df_raw["activity"] - y_hat
        else:
            self.df_raw["residual"] = self.df_raw["activity"] - np.mean(self.df_raw["activity"])

        t_test_result_1 = ttest_ind(
            x1=self.df_raw.query("treatment==1")["residual"],
            x2=self.df_raw.query("treatment==0")["residual"],
            alternative="larger",
        )  # test if peptide is an enhancer
        t_test_result_2 = ttest_ind(
            x1=self.df_raw.query("treatment==1")["residual"],
            x2=self.df_raw.query("treatment==0")["residual"],
            alternative="smaller",
        )  # test if peptide is a suppressor

        return AncovaResult(self.sample_id, significant_ph, round(t_test_result_1[1], 4), round(t_test_result_2[1], 4))
