#!/usr/bin/env python3

"""
This tool conducts a ANCOVA analysis for enzyme activity data. If no significant covariation detected between pH and activity, a t-test will be performed for verifying if the peptide has significant impact on activity. Otherwise a linear regression will be executed using pH as independent variable, and the residues will be used in the t-test.
"""

from __future__ import annotations

from pathlib import Path
from src.common import CovarProcessor
import argparse
import logging
import os


logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(filename)s %(message)s")
logging.getLogger().setLevel(logging.INFO)
logger = logging.getLogger(__file__)

def parse_args():
    """Parse any command line arguments.

    :return: Parsed arguments
    :rtype: argparse.Namespace
    """

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-id",
        "--sample-id",
        help="tested sample id ",
        type=str,
    )
    parser.add_argument(
        "-i",
        "--input-file",
        help="tsv file containing 3 columns: activity, ph, and treatment",
        type=Path,
    )
    parser.add_argument(
        "-p",
        "--p-value",
        help="maximal allowed P value for ANCOVA analysis",
        type=float,
    )
    
    return parser.parse_args()
    
def main(args: argparse.Namespace):
    
    with CovarProcessor(sample_id = args.sample_id, input_path = args.input_file, max_covar_p_value = args.p_value) as samplecovar:
        print(f"Running analysis for: {samplecovar.sample_id}")
        return samplecovar.t_test()
    
    
if __name__ == "__main__":
    args = parse_args()
    os.chdir("/home/jovyan/Documents/jupyter/pep-detective/")
    main(args)