PAIR2: efficient and powerful software to detect polymorphic Alu elements in high-throughput sequencing data. unpublished yet. Please contact me for usage

####
##############################
#### start from these files 
config.dk: include most of the paths of input/output files, parameters to use in the program. 
           used by Yu in Denmark, to test data mapped to hg19.
                  (@C: means mentioned parameters can be found in this config file)
qsub_dk.py: call programs or write commands to submit to cluster
	    most program has input like this:
            debug/program_to_call config.dk (option) (pn_idx)
	    eg: "debug/alu_insert config.dk alu_mate_flag 0" call alu_insert with 
                option "alu_mate_flag" for the 0th PN in the (@C)file_pn

config.decode: used by Yu in Iceland, where the old files sitting 

##################################
#### alu deletion 
1. build_dist.cpp
   /** create probability density function of insert length for each PN, stratified by reading group
       */
   results written into: (@C)file_dist_prefix 

2. alu_delete.cpp
   /** reads are classified as mid reads (cross alu region, indication of no deletion), clip reads (clip mapped at alu boundary, indication of 
       deletion) and unknown reads (can't know for sure if it is mid or clip read, use insert length instead) .   
       */
   opt = write_tmps_pn: 
       write temporary results: (@C)file_alu_delete0
       generate *log1: high coverage region
       generate *tmp1: counts of 3 types of read, and insert length for unknow reads. 
       generate *tmp2: calculte genotype prob based on *tmp1 files. 
                       header l0, l1, l2: genotype prob based on mid/clip reads
                       header m0, m1, m2: genotype prob based on mid/clip/unknow reads

   opt = write_vcf_pns:
       write final results: (@C)file_alu_delete0
       mutation that is private (only 1 individual has) is ignore. 
       generate *vcf and *pos files: as before, "lr" used mid/clip reads,  "mr" used mid/clip/unknow reads
       filter good quality mutations using *pos file.
       
   opt = debug:
       debug purpose
       
##################################
#### alu insertion (alu_insert.cpp)
1. alu_insert.cpp  (5 steps )
   opt = "write_tmps_pn": 
       write results: (@C)file_alu_insert0
       for each pn, scan the genome and find potention loci for alu insertion, this generate a lot of temporary files. 
       *tmp1 files can be ignored after this step
       *tmp1, reads mapped at dif chrn or far away
       *tmp1.chr1, result of *tmp1, where one end is mapped to alu 
       *tmp3, left_alu_cnt + right_alu_cnt >= 4 is considered. 
       *tmp3st, subset of *tmp3, with pos in repetitive regions are excluded 
       
   opt = "combine_pos_pns" :
        write results: (@C)file_alu_insert1
        find common insertion positions of multiple individuals, also filter out insertions outside pre_defined freq range.
        eg:  only insertions with freq in [minfreq, maxfreq] are considered, I use [0.02, 1] as in config.dk file

        Note: some PNs with extreme large amount of reads (that might be alu mate). (almost 10 times of other PNs)
              eg: MXQBTEO, HYKJEJS, EIVHBDQ, KSRDJGR
              EITHER:  > 4G memory to fun; OR: ignore them
              I chose to ignore them. For each chromosome, about 15% individuals with extremely high potential insert positions are ignored. 

	
   opt = "clipreads_pn":
     for each pn, write the broken reads at these insertion positions to files at (@C)file_clip_reads/chr1/pn_minfreq_maxfreq.
     a read is considered as broken reads only if cigar string shows soft clip, AND the score after realignment is high.

   opt = "clipreads_pns":
     write files at  (@C)file_clip_reads/chr1_pos/regionBegin_regionEnd
     all broken reads with insertion between regionBegin and regionEnd are written into corresponding files. 
     exact position is written at (@C)file_alu_insert1/chr1.clip_pn


   opt = "write_tmp0_pn":  	 
       	 see eg.sh

   opt = "clipPos_pns chr1":
         eg.sh     

   opt = "write_tmp2_pn":  	 
         eg.sh

   opt = "fixed_vcf_pns":  	 
     combine results and write big vcf 