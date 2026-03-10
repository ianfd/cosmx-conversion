from wf.conversion import cosmx_convert_with_stats_gen

from latch.resources.launch_plan import LaunchPlan
from latch.resources.workflow import workflow
from latch.types.directory import LatchDir, LatchOutputDir
from latch.types.file import LatchFile
from latch.types.metadata import LatchAuthor, LatchMetadata, LatchParameter


metadata = LatchMetadata(
    display_name="CosMx Conversion and Statistics",
    author=LatchAuthor(
        name="Ian",
    ),
    parameters={
        "sample_tar_gz": LatchParameter(
            display_name="Compressed CosMx Sample File",
            description="CosMx compressed flat file export (.tar.gz or similar).",
            batch_table_column=True,
        ),
        "sample_name": LatchParameter(
            display_name="Sample Name",
            description="Name used for output file naming.",
            batch_table_column=True,
        ),
        "output_dir": LatchParameter(
            display_name="Output Directory",
            description="Latch path for H5AD and statistics output.",
            batch_table_column=True,
        ),
    },
)


@workflow(metadata)
def cosmx_convert(
    sample_tar_gz: LatchFile,
    sample_name: str,
    output_dir: LatchDir = LatchDir("latch://40726.account/cosmx-test/out-dir"),
) -> LatchOutputDir:
    """
    ## CosMx Conversion + Statistics Generation

    Converts raw CosMx flat file exports into H5AD format
    and generates pre-QC statistics for inspection in Latch Plots.

    ### Steps
    1. Extract and load CosMx data with squidpy
    2. Compute QC metrics (total counts, genes detected, negative probes, protein markers)
    3. Generate per-cell, per-FOV, protein, and summary statistics CSVs
    4. Save H5AD with raw counts preserved

    ### Outputs
    - `{sample_name}.h5ad` — raw AnnData with QC metrics in obs
    - `stats/` — pre-filter statistics CSVs for visualization in Plots
    """
    return cosmx_convert_with_stats_gen(
        sample_tar_gz=sample_tar_gz,
        sample_name=sample_name,
        output_dir=output_dir,
    )



LaunchPlan(
    cosmx_convert,
    "Test Data",
    {
        "sample_tar_gz": LatchFile("latch://40726.account/cosmx-test/GSE282193_Slide1.tar.gz"),
        "sample_name": "GSE282193_Slide1",
        "output_dir": LatchDir("latch://40726.account/cosmx-test/out-dir/conv/GSE282193_Slide1"),
    }
)