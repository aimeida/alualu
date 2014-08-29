import sys, os

def print2(pn_all, fn_path, path1, bin_path, time_sec):
    if not os.path.exists(fn_path):
        os.mkdir(fn_path)
    pi = 0
    for pn in pn_all:
        with open(fn_path + str(pi), 'w') as fout:
            print >>fout, "#!/bin/sh "
            print >>fout, "#SBATCH -N 1"  ## nodes 
            print >>fout, "#SBATCH -c 1" ## ppn
            if time_sec < 60:
                print >>fout, "#SBATCH -p express"
                print >>fout, "#SBATCH -t 0-0:%d:0"%time_sec
            else:
                print >>fout, "#SBATCH -p normal"
                print >>fout, "#SBATCH -t 0-8:0:0" 

            print >>fout, "#SBATCH --mem=4g" 
            print >>fout, "#SBATCH --job-name %d_%s" % (pi, pn)
            print >>fout, "#SBATCH -o %d_%s.o" % (pi, pn)
            print >>fout, "#SBATCH -e %d_%s.e" % (pi, pn)
            print >>fout, "echo `date`"
            print >>fout, '%(bin_path)salu_delete ../config.dk write_tmps_pn %(pi)d'%locals()
            #print >>fout, '%(bin_path)salu_insert ../config.dk clipReads_pn %(pi)d'%locals()
            #print >>fout, '%(bin_path)salu_insert ../config.dk write_tmp0_pn %(pi)d'%locals()
            #print >>fout, '%(bin_path)salu_insert ../config.dk write_tmp2_pn %(pi)d'%locals()
            print >>fout, "echo `date`"

        pi += 1


def print3(fn_path, path1, bin_path):
    if not os.path.exists(fn_path):
        os.mkdir(fn_path)
    pi = 0
    for _chrn in range(1,23) + ['X', 'Y']:
        chrn = "chr" + str(_chrn)
        with open(fn_path + chrn, 'w') as fout:
            print >>fout, "#!/bin/sh "
            print >>fout, "#SBATCH -p express"
            print >>fout, "#SBATCH -N 1"  ## nodes 
            print >>fout, "#SBATCH -c 1" ## ppn
            print >>fout, "#SBATCH -t 0-0:59:0"
            print >>fout, "#SBATCH --mem=4g" 
            print >>fout, "#SBATCH --job-name %s" % (chrn)
            print >>fout, "#SBATCH -o %s.o" % chrn
            print >>fout, "#SBATCH -e %s.e" % chrn
            #print >>fout, '%(bin_path)salu_insert ../config.dk clipReads_pns %(chrn)s'%locals()
            print >>fout, '%(bin_path)salu_insert ../config.dk clipPos_pns %(chrn)s'%locals()
        pi += 1


def ins_del(fn, bin_path, opt):
    with open(fn,'w') as fout:
        print >> fout, '%(bin_path)s%(opt)s/alu_delete config.dk write_vcf_pns'%locals()
        print >> fout, '%(bin_path)s%(opt)s/alu_insert config.dk preprocess'%locals()
        print >> fout, '%(bin_path)s%(opt)s/alu_insert config.dk combine_pos_pns'%locals()
        print >> fout, '%(bin_path)s%(opt)s/alu_insert config.dk fixed_vcf_pns'%locals()
    print 'written into', fn


if __name__ == '__main__':
    path1 = '/home/qianyuxx/faststorage/Alu/'
    pn_all = map(lambda x:x.strip(), file('/home/qianyuxx/faststorage/AluDK/inputs/PN_all').readlines())
    print2(pn_all, path1+"q_ad/", path1, path1 + 'opt3/', 159)
    #print2(pn_all, path1+"q_ad/", path1, path1 + 'opt3/', 59)

    pn_used = map(lambda x:x.strip(), file('/home/qianyuxx/faststorage/AluDK/outputs/insert_alu1/pn_used').readlines())
    #print2(pn_used, path1+"q_ai/", path1, path1 + 'opt3/', 59)
    #print2(pn_used, path1+"q_ai/", path1, path1 + 'debug/', 159)
    
    ## run by chr
    #print3(path1+"q_ai/", path1, path1 + 'opt3/')
    #print3(path1+"q_ai/", path1, path1 + 'debug/')

    #ins_del("ins_del.sh", path1, "opt3")
