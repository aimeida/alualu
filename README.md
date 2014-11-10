PAIR2: efficient and powerful software to detect polymorphic Alu elements in high-throughput sequencing data. 
NB:   unpublished yet. Please contact me for usage
configFile: indicate the paths of input/output files and parameters to use in the program. 
@p: denote parameter 'p' given in configFile

###################################
##### run alu_deletion (make opt3/alu_delete)
     (1) run "alu_delete configFile preprocess"
     	 sort the alu position file, make folders for outputs
 
     (2) run "alu_delete configFile write_tmps_pn X" 
         X is index of individual according to the @file_pn. eg. for 2000 individuals, we can have 2000 jobs running in parallel, X ranging from 0 to 1999. 
	 30X sequencing depth takes about 5hrs per individual (very rarely it may exceed 12 hrs)  

     (3) run "alu_delete configFile write_vcf_pns"
         make big vcf files     

###################################
##### run alu_insertion (make opt3/alu_insert)
      step(2) takes a few hours to run, all the other steps takes < 60 min, < 4G memory is enough for 99.9% samples.

     (1)  alu_insert configFile preprocess  

     (2)  run "alu_insert configFile write_tmps_pn X", X is index of an individual.

     (3)  make a file, @file_pn_used 
	  include a subset of individuals from file @file_pn, individuals failed to run step (2) or with large amount of discordant reads are excluded.
	  eg. 2000 individuals, up to 10% can be removed because it takes too much time to handle their discordant reads.

     (4) run "alu_insert configFile combine_pos_pns"
          which generate files named @file_alu_insert1/insert_pos.chr*  

     (5) run "alu_insert configFile clipReads_pn X" for each pn
         write positions of clip reads at @file_alu_insert1/clip/chr*/pn
	 X is index of individual, if its corresponding name does not exist in @file_pn_used, the program will quit

     (6) run "alu_insert configFile clipReads_pns chr*" for each chromosome (chr1 - chrX)
     	 write insert_alu1/clip/chr*_pos/*
         which generate files at  insert_alu1/clip/chr*.clip_pn and  insert_alu1/clip/chr*.clip_region
	
     (7) run "alu_insert configFile write_tmp0_pn X" for each pn
         which writes alu and clip reads for each individual, at insert_alu1/cons/chr*/*

     (8) run "alu_insert configFile clipPos_pns chrX" for each chromosome (chr1 - chrX)
         write exact insertion breakpoints at  insert_alu1/cons/chr*_clip_pass       

     (9) run "alu_insert configFile write_tmp2_pn X" for each pn
         genotype calling      

     (10) run "alu_insert configFile fixed_vcf_pns" 
         final vcf file
