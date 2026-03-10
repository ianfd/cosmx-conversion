from wf.__init__ import *

cosmx_convert_with_stats_gen(
        expr_mat= LatchFile("latch://40726.account/cosmx-test/GSE282193_Slide1.tar.gz"),
        sample_name= "GSE282193_Slide1",
        output_dir= LatchDir("latch://40726.account/cosmx-testf/out-dir/conv/GSE282193_Slide1"),
)