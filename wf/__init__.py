from typing import Dict, List

from wf.conversion import cosmx_convert_with_stats_gen

from latch.resources.launch_plan import LaunchPlan
from latch.resources.workflow import workflow
from latch.types.directory import LatchDir, LatchOutputDir
from latch.types.file import LatchFile
from latch.types.metadata import LatchAuthor, LatchMetadata, LatchParameter
from latch.resources.tasks import small_task, medium_task
from latch.resources.map_tasks import map_task 
import extra.loading_utils


metadata = LatchMetadata(
    display_name="CosMx Conversion and Statistics for multiple samples",
    author=LatchAuthor(
        name="Ian",
    ),
    parameters={
        "sample_files": LatchParameter(
            display_name="List of Sample Files ",
            description="Dictionary of Sample Files",
            batch_table_column=True,
        ),
        "output_dir_base": LatchParameter(
            display_name="Output Directory",
            description="Latch path for H5AD and statistics output.",
            batch_table_column=True,
        ),
    },
)


# @workflow(metadata)
# def cosmx_convert(
#     sample_tar_gz: LatchFile,
#     sample_name: str,
#     output_dir: LatchDir = LatchDir("latch://40726.account/cosmx-test/out-dir"),
# ) -> LatchOutputDir:
#     """
#     ## CosMx Conversion + Statistics Generation

#     Converts raw CosMx flat file exports into H5AD format
#     and generates pre-QC statistics for inspection in Latch Plots.

#     ### Steps
#     1. Extract and load CosMx data with squidpy
#     2. Compute QC metrics (total counts, genes detected, negative probes, protein markers)
#     3. Generate per-cell, per-FOV, protein, and summary statistics CSVs
#     4. Save H5AD with raw counts preserved

#     ### Outputs
#     - `{sample_name}.h5ad` — raw AnnData with QC metrics in obs
#     - `stats/` — pre-filter statistics CSVs for visualization in Plots
#     """
#     return cosmx_convert_with_stats_gen(
#         sample_tar_gz=sample_tar_gz,
#         sample_name=sample_name,
#         output_dir=output_dir,
#     )


@small_task
def prep_args_for_multi(sample_files: List[str],
    sample_names: List[str], base_dir: LatchDir) -> List[extra.loading_utils.ConvertInput]:
    
    return [extra.loading_utils.ConvertInput(file=LatchFile(v), output_dir=base_dir, sample_name=k) for k,v in zip(sample_names, sample_files)]

@workflow(metadata)
def cosmx_convert_multi(
    sample_files: List[str],
    sample_names: List[str],
    output_dir_base: LatchDir = LatchDir("latch://40726.account/cosmx-test/out-dir"),
) -> List[LatchOutputDir]:

    prepped_args = prep_args_for_multi(sample_files=sample_files, sample_names=sample_names, base_dir=output_dir_base)

    mapped_out = map_task(cosmx_convert_with_stats_gen)(input=prepped_args)
    return mapped_out


LaunchPlan(
    cosmx_convert_multi,
    "Test Data",
    {
        "sample_names": ["GSE282193_Slide1"],
        "sample_files": ["latch://40726.account/cosmx-test/GSE282193_Slide1.tar.gz"],
        "output_dir_base": LatchDir("latch://40726.account/cosmx-testf/out-dir/conv/"),
    }
)