
make OD=opt3 opt3/alu_insert2
#opt3/alu_insert2 config.dk cons_reads_build chr1 ## running now 


# opt3/alu_delete config.dk write_tmps_pn 0
# opt3/alu_delete config.dk write_vcf_pns

# opt3/alu_insert config.dk write_tmps_pn 0
# opt3/alu_insert config.dk combine_pos_pns
# opt3/alu_insert config.dk clipReads_by_pn 0
# opt3/alu_insert config.dk clipReads_pos_pns
# opt3/alu_insert config.dk fixed_delete0_pn 0
# opt3/alu_insert config.dk fixed_vcf_pns 

make OD=opt3 opt3/alu_insert
#debug/alu_insert config.dk clipReads_pos_pns

