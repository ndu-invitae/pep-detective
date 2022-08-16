#!/usr/bin/env python3

import logging
from pathlib import Path
from typing import Union

import pandas as pd

from pep_detective.src.conn import get_db, init_conn


def get_all_db(dbs):
    return pd.concat([get_db(init_conn(db)) for db in dbs])


get_all_db(snakemake.input).to_csv(snakemake.output[0], sep="\t", index=False)
