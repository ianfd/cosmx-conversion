from wf.__init__ import *

cosmx_convert_with_stats_gen(
        sample_names= ["GSE282193_Slide1"],
        sample_files= ["latch://40726.account/cosmx-test/GSE282193_Slide1.tar.gz"],
        output_dir_base = LatchDir("latch://40726.account/cosmx-testf/out-dir/conv/"),
)