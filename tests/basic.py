from wf.__init__ import *

cosmx_convert_with_stats_gen(
        samples = {"GSE282193_Slide1": LatchFile("latch://40726.account/cosmx-test/GSE282193_Slide1.tar.gz")},
        output_dir_base = LatchDir("latch://40726.account/cosmx-testf/out-dir/conv/"),
)