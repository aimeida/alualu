opt="opt3"
${opt3}/alu_delete config.dk write_tmps_pn 0
${opt3}/alu_delete config.dk write_vcf_pns

${opt3}/alu_insert config.dk write_tmps_pn 0
${opt3}/alu_insert config.dk combine_pos_pns     ## need to provide file_pn_used
${opt3}/alu_insert config.dk clipReads_pn 0      ## start from here, remove bad reads 
${opt3}/alu_insert config.dk clipReads_pns chr1 
${opt3}/alu_insert config.dk write_tmp0_pn 0
${opt3}/alu_insert config.dk clipPos_pns chr1  ## find exact clipPos
${opt3}/alu_insert config.dk write_tmp2_pn 0  ## first count clip and alu reads, then genotype calling 
${opt3}/alu_insert config.dk fixed_vcf_pns

