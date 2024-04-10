r, pt, st, sz = [800, .01, .85, 1]
params = f"{r}__{pt}__{st}__{sz}__mustache"
path = f"final_loops"
file='../cools/DOWNSAMPLED_merged_all_800.cool'
cmd = f"mkdir -p final_loops; mkdir -p mustache_output; mustache -f {file} -o ./{path}/out_merged_DOWNSAMPLE.tsv -p 8 -pt {pt} -st {st} -sz {sz} -r {r}"

import subprocess
subprocess.run(cmd, shell=True)
