######## MAKE MERGED
pref=../cools/
files=(nc1-8_1_1+1_2+2_1+2_2.base_resolution=100.mcool nc14_1_1+1_2+2_1+2_2.base_resolution=100.mcool s10-12_1_1+1_2+2_1+2_2.base_resolution=100.mcool)

for file in ${files[@]}
do
        echo doing ${file}
        cooltools random-sample -p 16 -c 287727300 --exact  ${pref}/${file}::/resolutions/800  ./DOWNSAMPLED_${file}_800
        cooler balance ./DOWNSAMPLED_${file}_800
        echo done with ${file}
done

cooler merge ../cools/DOWNSAMPLED_merged_all_800.cool DOWNSAMPLED_nc1-8_1_1+1_2+2_1+2_2.base_resolution=100.mcool_800 DOWNSAMPLED_nc14_1_1+1_2+2_1+2_2.base_resolution=100.mcool_800 DOWNSAMPLED_s10-12_1_1+1_2+2_1+2_2.base_resolution=100.mcool_800
cooler balance ../cools/DOWNSAMPLED_merged_all_800.cool





######### COMPUTE EXPECTED
files=(nc14_1_1.base_resolution=100.mcool nc1-8_1_2.base_resolution=100.mcool nc14_1_2.base_resolution=100.mcool nc1-8_2_1.base_resolution=100.mcool  nc14_1_3.base_resolution=100.mcool nc1-8_2_2.base_resolution=100.mcool nc14_2_1.base_resolution=100.mcool nc14_2_2.base_resolution=100.mcool s10-12_1_1.base_resolution=100.mcool  s10-12_1_2.base_resolution=100.mcool  s10-12_2_1.base_resolution=100.mcool  nc1-8_1_1.base_resolution=100.mcool s10-12_2_2.base_resolution=100.mcool nc14_FED_1_1.base_resolution=100.mcool nc14_FED_1_2.base_resolution=100.mcool  nc14_FED_rerun_1_1.base_resolution=100.mcool nc14_FED_rerun_1_2.base_resolution=100.mcool nc12_mitotic_FED_rerun_1_1.base_resolution=100.mcool nc12_mitotic_FED_rerun_1_2.base_resolution=100.mcool nc12_mitotic_FED_1_1.base_resolution=100.mcool nc12_mitotic_FED_1_2.base_resolution=100.mcool)
for file in ${files[@]}
do
       echo doing ${file}
       cooltools compute-expected ${pref}/${file}::/resolutions/800 -p 16 -o expected/expected_${file}_800.tsv
       echo done with ${file}
done

cooltools compute-expected ${pref}/DOWNSAMPLED_merged_all_800.cool -p 16 -o  expected/expected_DOWNSAMPLED_merged_all_800.cool.tsv








