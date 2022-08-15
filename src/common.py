from __future__ import annotations

import pandas as pd
import numpy as np
from pingouin import ancova
from statsmodels.formula.api import ols
from statsmodels.stats.weightstats import ttest_ind
from pathlib import Path
from dataclasses import dataclass


@dataclass
class AncovaResult:
    """Dataclass for storing ANCOVA analysis results."""
    sample_id:str
    ph_covar:bool
    p_enhancer:float
    p_suppressor:float    
    
class CovarProcessor:
    """
    Main class for processing activity data.
    

    input_path: Input tsv file that contains one dependent variable (dv) column, one between-subjects factor column with two levels (bf), and one covariate (covar) column
    max_covar_p_value: maximal allowed P value for rejecting no significant correlation hypothesis (default 0.05)
    """
    
    def __init__(
        self,
        sample_id:str,
        input_path:Path,
        max_covar_p_value:float = 0.05,
        ):
        
        self.sample_id = sample_id
        self.input_path = input_path
        self.max_covar_p_value = max_covar_p_value

        self._open_tsv()
        
        def _open_tsv(self):
            """load in tsv file"""
            
            self.df_raw = pd.read_table(self.input_path)

        
        @property    
        def ancova_analysis(self) -> bool:
            """Run ANCOVA and determine if pH is a significant covariate """
            
            covar_stats = ancova(data = self.df_raw, dv = 'activity',between="treatment", covar="ph")
            self.significant_ph = (covar_stats.query('Source == "ph"')['p-unc'] < self.max_covar_p_value)
            
        def t_test(self) -> AncovaResult:
            """Execute t test on activities with/ without peptide treatment and return """
            
            if self.significant_ph: # perform linear fitting if pH is a significant covariate otherwise use mean for residue calculaition
                lm = ols("activity ~ ph", self.df_raw).fit()
                y_hat = lm.predict(self.df_raw['ph'])
                self.df_raw["residual"] = self.df_raw['activity'] - y_hat
            else:
                self.df_raw["residual"] = self.df_raw['activity'] - np.mean(self.df_raw['activity'])
                
            t_test_result_1 = ttest_ind(x1=self.df_raw.query('treatment==1')['residual'], x2=self.df_raw.query('treatment==0')['residual'], alternative = "larger") # test if peptide is an enhancer
            t_test_result_2 = ttest_ind(x1=self.df_raw.query('treatment==1')['residual'], x2=self.df_raw.query('treatment==0')['residual'], alternative = "smaller") # test if peptide is a suppressor
            
            return AncovaResult(self.sample_id, self.significant_ph, t_test_result_1[1], t_test_result_2[1])
            
    
    