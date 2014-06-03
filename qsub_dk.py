import sys, os

def print1(pn_all, fn_path, path1):
    pi = 0
    for pn in pn_all:
        with open(fn_path + str(pi), 'w') as fout:
            print >>fout, "#!/bin/sh "
            print >>fout, "#PBS -q normal"
            print >>fout, "#PBS -l nodes=1:ppn=3"
            print >>fout, "#PBS -N %s" % pn
            print >>fout, "cd %s" % path1
            #print >>fout, 'debug/build_dist config.dk %(pi)d'%locals() 
            print >>fout, 'debug/alu_delete config.dk write_tmp1 %(pi)d'%locals() 
            print >>fout, 'debug/alu_delete config.dk write_tmp2 %(pi)d'%locals() 

            #print >>fout, 'debug/alu_insert config.dk alu_mate_flag %(pi)d'%locals()
            #print >>fout, 'debug/alu_insert config.dk combine_pn_pos %(pi)d'%locals()

            #print >>fout, 'debug/insert_pos config.dk 1 %(pi)d'%locals()
            print >>fout, 'debug/ins_del config.dk write_tmp1 %(pi)d'%locals() 
        pi += 1

if __name__ == '__main__':
    path1 = os.getcwd() + '/'
    pn_all = map(lambda x:x.strip(), file(path1 + 'inputs/PN_all').readlines())
    #print1(pn_all, path1+"q_dist/", path1)
    #print1(pn_all, path1+"q_ai/", path1)
    print1(pn_all, path1+"q_ad/", path1)
