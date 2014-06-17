import sys, os

def print1(pn_all, fn_path, path1, bin_path, fastq = True):
    if not os.path.exists(fn_path):
        os.mkdir(fn_path)
    pi = 0
    for pn in pn_all:
        with open(fn_path + str(pi), 'w') as fout:
            print >>fout, "#!/bin/sh "
            print >>fout, "#PBS -q normal"
            print >>fout, "#PBS -l nodes=1:ppn=1"
            if fastq:
                print >>fout, "#PBS -l walltime=0:59:0"
            print >>fout, "#PBS -N %d_%s" % (pi, pn)
            print >>fout, "cd %s" % path1
            #print >>fout, '%(bin_path)sbuild_dist config.dk build_dist %(pi)d'%locals() 
            #print >>fout, '%(bin_path)salu_delete config.dk write_tmps %(pi)d'%locals() 
            print >>fout, '%(bin_path)salu_insert config.dk write_tmps %(pi)d'%locals()
            #print >>fout, '%(bin_path)salu_insert config.dk combine_pos %(pi)d'%locals()            
            #print >>fout, '%(bin_path)sinsert_pos config.dk clipReads_by_pn %(pi)d'%locals()
            #print >>fout, '%(bin_path)sinsert_pos config.dk cons_reads_pn %(pi)d'%locals()
        pi += 1

if __name__ == '__main__':
    path1 = os.getcwd() + '/'
    pn_all = map(lambda x:x.strip(), file('/home/qianyuxx/faststorage/AluDK/inputs/PN_all').readlines())
    #print1(pn_all, path1+"q_dist/", path1)
    #print1(pn_all, path1+"q_ai/", path1, path1 + 'debug/', True)
    print1(pn_all, path1+"q_ai/", path1, path1 + 'opt3/', True)
    #print1(pn_all, path1+"q_ai/", path1, path1 + 'opt3/', False)
    #print1(pn_all, path1+"q_ad/", path1, path1 + 'opt3/', False)
