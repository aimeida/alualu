import sys, os

def print1(pn_all, fn_path, path1, bin_path, fast_queue = True):
    if not os.path.exists(fn_path):
        os.mkdir(fn_path)
    pi = 0
    for pn in pn_all:
        with open(fn_path + str(pi), 'w') as fout:
            print >>fout, "#!/bin/sh "
            print >>fout, "#PBS -q normal"
            print >>fout, "#PBS -l nodes=1:ppn=1"
            if fast_queue:
                print >>fout, "#PBS -l walltime=0:59:0"
            print >>fout, "#PBS -N %d_%s" % (pi, pn)
            print >>fout, "cd %s" % path1
            #print >>fout, '%(bin_path)sbuild_dist config.dk build_dist %(pi)d'%locals() 
            #print >>fout, '%(bin_path)salu_delete config.dk write_tmps_pn %(pi)d'%locals() 
        pi += 1

def print2(pn_all, fn_path, path1, bin_path, fast_queue = True):
    if not os.path.exists(fn_path):
        os.mkdir(fn_path)
    pi = 0
    for pn in pn_all:
        with open(fn_path + str(pi), 'w') as fout:
            print >>fout, "#!/bin/sh "
            print >>fout, "#SBATCH -p normal"
            print >>fout, "#SBATCH -N 1"  ## nodes 
            print >>fout, "#SBATCH -c 1" ## ppn
            if fast_queue:
                print >>fout, "#SBATCH -t 0-0:59:0"
            print >>fout, "#SBATCH --mem=4g" 
            print >>fout, "#SBATCH --job-name %d_%s" % (pi, pn)
            print >>fout, "#SBATCH -o %d_%s.o" % (pi, pn)
            print >>fout, "#SBATCH -e %d_%s.e" % (pi, pn)
            print >>fout, '%(bin_path)salu_delete ../config.dk write_tmps_pn %(pi)d'%locals()
            #print >>fout, '%(bin_path)salu_insert ../config.dk clipReads_pn %(pi)d'%locals()
            #print >>fout, '%(bin_path)salu_insert ../config.dk write_tmp0_pn %(pi)d'%locals()
            #print >>fout, '%(bin_path)salu_insert ../config.dk write_tmp2_pn %(pi)d'%locals()
        pi += 1


def print3(fn_path, path1, bin_path):
    if not os.path.exists(fn_path):
        os.mkdir(fn_path)
    pi = 0
    for _chrn in range(1,23) + ['X', 'Y']:
        chrn = "chr" + str(_chrn)
        with open(fn_path + chrn, 'w') as fout:
            print >>fout, "#!/bin/sh "
            print >>fout, "#SBATCH -p normal"
            print >>fout, "#SBATCH -N 1"  ## nodes 
            print >>fout, "#SBATCH -c 1" ## ppn
            print >>fout, "#SBATCH -t 0-0:59:0"
            print >>fout, "#SBATCH --mem=4g" 
            print >>fout, "#SBATCH --job-name %s" % (chrn)
            print >>fout, "#SBATCH -o %s.o" % chrn
            print >>fout, "#SBATCH -e %s.e" % chrn
            #print >>fout, '%(bin_path)salu_insert ../config.dk clipReads_pns %(chrn)s'%locals()
            #print >>fout, '%(bin_path)salu_insert ../config.dk clipPos_pns %(chrn)s'%locals()
            #print >>fout, '%(bin_path)salu_insert2 ../config.dk consReads_build %(chrn)s'%locals()
            #print >>fout, '%(bin_path)salu_insert2 ../config.dk consReads_chr %(chrn)s'%locals()
        pi += 1

if __name__ == '__main__':
    path1 = os.getcwd() + '/'
    pn_all = map(lambda x:x.strip(), file('/home/qianyuxx/faststorage/AluDK/inputs/PN_all').readlines())
    #print1(pn_all, path1+"q_dist/", path1)

    pn_all = map(lambda x:x.strip(), file('/home/qianyuxx/faststorage/AluDK/outputs/insert_alu1/pn_used').readlines())
    print2(pn_all, path1+"q_ad/", path1, path1 + 'opt3/', True)
    #print2(pn_all, path1+"q_ai/", path1, path1 + 'opt3/', True)
    
    ## run by chr
    #print3(path1+"q_ai/", path1, path1 + 'opt3/')
    
