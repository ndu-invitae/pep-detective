import json
from pathlib import Path

import numpy as np
import pandas as pd

np.random.seed(666)


def create_fake_data(num_sample: int = 50, num_rep: int = 5):
    """create fake tsv files and a fake json file that minics LIMS query result; each tsv file is coresponding to one sample."""

    def _activity(x, a, b, e):
        return a * x + b + e

    def _retrive_sample():
        for sample in df_raw.groupby("sample"):
            yield sample

    temp_path = Path("tmp")
    temp_path.mkdir(parents=True, exist_ok=True)

    ph = np.random.normal(loc=0, scale=0.2, size=num_sample * num_rep * 2) + 7
    slope = np.random.normal(loc=0.5, scale=0.2, size=num_sample)
    factor = np.array([1] * num_rep + [0] * num_rep)
    intercept = np.transpose(
        np.random.gamma(shape=1.2, scale=5, size=num_sample * num_rep * 2).reshape(num_sample, num_rep * 2) * factor
    )
    error = np.random.normal(loc=0, scale=0.5, size=num_sample * num_rep * 2).reshape(num_rep * 2, num_sample) + 3

    df_raw = (
        pd.DataFrame(
            _activity(x=ph.reshape(num_rep * 2, num_sample), a=slope, b=intercept, e=error),
            index=factor,
            columns=[f"smaple{i+1}" for i in range(num_sample)],
        )
        .stack()
        .reset_index()
    )
    df_raw.columns = ["treatment", "sample", "activity"]
    df_raw["ph"] = ph
    df_raw = df_raw.set_index(["sample"])

    dict_json = {}
    for sample in _retrive_sample():
        path = temp_path / f"{sample[0]}.tsv"
        dict_json.update({sample[0]: str(path.resolve())})
        sample[1][["treatment", "ph", "activity"]].to_csv(path, sep="\t", index=False)

    with open(temp_path / "input_files.json", "w") as fo:
        json.dump(dict_json, fo)
