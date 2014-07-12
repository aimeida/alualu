opt="opt3"
${opt3}/alu_delete config.dk write_tmps_pn 0
${opt3}/alu_delete config.dk write_vcf_pns

${opt3}/alu_insert config.dk write_tmps_pn 0
${opt3}/alu_insert config.dk combine_pos_pns     ## need to provide file_pn_used
${opt3}/alu_insert config.dk clipReads_pn 0      ## can also start from here, remove bad reads 
${opt3}/alu_insert config.dk clipReads_pns chr1  ## start from here  
${opt3}/alu_insert config.dk write_tmp0_pn 0
${opt3}/alu_insert config.dk write_tmp1 chr1  ## clip and alu reads
${opt3}/alu_insert config.dk write_tmp2_pn 0  ## clip and alu reads
${opt3}/alu_insert config.dk fixed_vcf_pns

