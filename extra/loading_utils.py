import os
import re
from pathlib import Path
from typing import TYPE_CHECKING

from pathlib import Path
import tarfile
import gzip
import shutil

from latch.types.directory import LatchDir, LatchOutputDir
from latch.types.file import LatchFile
import io

import pandas as pd
from anndata import AnnData
from scipy.sparse import csr_matrix

from dataclasses import dataclass


@dataclass
class ConvertInput:
    sample_name: str
    file: LatchFile
    output_dir: LatchDir

if TYPE_CHECKING:
    from dask.dataframe import DataFrame as DaskDataFrame


class CosmxKeys:
    COUNTS_SUFFIX = "exprMat_file.csv"
    TRANSCRIPTS_SUFFIX = "tx_file.csv"
    METADATA_SUFFIX = "metadata_file.csv"
    FOV_SUFFIX = "fov_positions_file.csv"
    FOV = "fov"
    CELL_ID = "cell_ID"
    INSTANCE_KEY = "cell"
    REGION_KEY = "fov_labels"
    X_GLOBAL_CELL = "CenterX_global_px"
    Y_GLOBAL_CELL = "CenterY_global_px"
    X_LOCAL_CELL = "CenterX_local_px"
    Y_LOCAL_CELL = "CenterY_local_px"
    X_LOCAL_TRANSCRIPT = "x_local_px"
    Y_LOCAL_TRANSCRIPT = "y_local_px"
    TARGET_OF_TRANSCRIPT = "target"


def cosmx_simple(
    path: str | Path,
    sample_name: str,
    dataset_id: str | None = None,
    *,
    counts_df: pd.DataFrame | None = None,
    obs_df: pd.DataFrame | None = None,
    fov_positions_df: pd.DataFrame | None = None,
) -> AnnData:
    import gc
    import numpy as np

    path = Path(path)

    fov_key = CosmxKeys.FOV
    cell_id_key = CosmxKeys.CELL_ID
    instance_key = CosmxKeys.INSTANCE_KEY

    if counts_df is None or obs_df is None:
        if dataset_id is None:
            counts_files = [f for f in os.listdir(path) if f.endswith(CosmxKeys.COUNTS_SUFFIX)]
            if len(counts_files) == 1:
                found = re.match(rf"(.*)_{CosmxKeys.COUNTS_SUFFIX}", counts_files[0])
                if found:
                    dataset_id = found.group(1)
        if dataset_id is None:
            raise ValueError("Could not infer `dataset_id` from the name of the counts file. Please specify it manually.")

        if counts_df is None:
            counts_df = pd.read_csv(
                path / f"{dataset_id}_{CosmxKeys.COUNTS_SUFFIX}",
                header=0, dtype=np.float32,
            )
        if obs_df is None:
            obs_df = pd.read_csv(
                path / f"{dataset_id}_{CosmxKeys.METADATA_SUFFIX}",
                header=0,
            )
        if fov_positions_df is None:
            fov_file = path / f"{dataset_id}_{CosmxKeys.FOV_SUFFIX}"
            if fov_file.exists():
                fov_positions_df = pd.read_csv(fov_file, header=0, index_col=fov_key)

    cell_map = obs_df[[cell_id_key, fov_key, instance_key]].drop_duplicates()
    counts_df = counts_df.merge(cell_map, on=[cell_id_key, fov_key], how="left")
    del cell_map
    gene_cols = [c for c in counts_df.columns if c not in [cell_id_key, fov_key, instance_key]]
    counts_df = counts_df.set_index(instance_key)[gene_cols]

    for col in obs_df.columns:
        s = obs_df[col]
        if isinstance(s.dtype, pd.CategoricalDtype):
            obs_df[col] = s.astype(str)
        elif pd.api.types.is_string_dtype(s):
            obs_df[col] = s.to_numpy(dtype=str, na_value="")
        elif hasattr(s.dtype, "numpy_dtype"):
            obs_df[col] = s.to_numpy(dtype=s.dtype.numpy_dtype)

    obs_df[fov_key] = pd.Categorical(obs_df[fov_key].astype(str))
    obs_df = obs_df.set_index(instance_key)
    obs_df.index = obs_df.index.astype(object)
    counts_df.index = counts_df.index.astype(object)

    common_index = obs_df.index.intersection(counts_df.index)

    counts_values = counts_df.loc[common_index, :].values
    gene_names = counts_df.columns.astype(object)
    del counts_df
    gc.collect()

    X = csr_matrix(counts_values)
    del counts_values
    gc.collect()

    obs_subset = obs_df.loc[common_index, :]
    del obs_df

    adata = AnnData(X, obs=obs_subset)
    del X, obs_subset
    adata.var_names = gene_names

    adata.obsm["spatial"] = adata.obs[["CenterX_local_px", "CenterY_local_px"]].values
    adata.obsm["spatial_fov"] = adata.obs[["CenterX_global_px", "CenterY_global_px"]].values
    adata.obs.drop(columns=["CenterX_local_px", "CenterY_local_px", "CenterX_global_px", "CenterY_global_px"], inplace=True)

    if fov_positions_df is not None:
        for fov in adata.obs[fov_key].cat.categories:
            adata.uns.setdefault("spatial", {})[fov] = {
                "images": {},
                "scalefactors": {"tissue_hires_scalef": 1, "spot_diameter_fullres": 1},
            }
        for fov, row in fov_positions_df.iterrows():
            if str(fov) in adata.uns.get("spatial", {}):
                adata.uns["spatial"][str(fov)]["metadata"] = row.to_dict()

    adata.obs['sample'] = sample_name

    return adata


def extract_files(input_file: str, work_dir: str | Path) -> str:
    sample_path = Path(input_file)   
    work_path = Path(work_dir)

    temp_dir = work_path / f"{sample_path.name.replace('.tar.gz', '')}_tmp"
    tar_extract_dir = temp_dir / "tar_contents"
    tar_extract_dir.mkdir(parents=True, exist_ok=True)

    with tarfile.open(sample_path, "r:gz") as tar:
        tar.extractall(path=tar_extract_dir)

    csv_gz_files = list(tar_extract_dir.rglob("*.csv.gz"))
    for gz_file in csv_gz_files:
        out_file = (temp_dir / gz_file.relative_to(tar_extract_dir)).with_suffix("")
        out_file.parent.mkdir(parents=True, exist_ok=True)

        with gzip.open(gz_file, "rb") as f_in, open(out_file, "wb") as f_out:
            shutil.copyfileobj(f_in, f_out)

    # Remove raw tar contents to free disk space
    shutil.rmtree(tar_extract_dir, ignore_errors=True)

    return str(temp_dir)

def prep_fov_file(sample_folder: str):
    target_file = ""

    for f in os.listdir(sample_folder):
        #print(f) debug
        if f.endswith("_fov_positions_file.csv"):
            target_file = f
    
    if target_file == "":
        return
    
    smpl_full_path = str(sample_folder) + "/" + target_file

    d = pd.read_csv(smpl_full_path)

    d = d.rename({"FOV": "fov"}, axis="columns")
    
    d.to_csv(smpl_full_path)


def read_full_sample(sample_name: str, data: str, work_dir: str | Path):
    temp_dir = extract_files(data, work_dir)
    prep_fov_file(temp_dir)

    spatial_data = cosmx_simple(
        path=temp_dir,
        sample_name=sample_name,
    )

    shutil.rmtree(temp_dir, ignore_errors=True)

    return spatial_data    