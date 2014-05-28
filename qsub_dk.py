import sys, os

### qx -v --nodes=5 -i inputs/bamFile/* -i qsub0.txt dispatch -r qsub0.txt > tmp3.sh
### qx --nodes=5 -i qsub0.txt dispatch -r qsub0.txt

def print1(pn_all, fn_path, path1):
    pi = 0
    for pn in pn_all:
        with open(fn_path + str(pi), 'w') as fout:
            print >>fout, "#!/bin/sh "
            print >>fout, "#PBS -q normal"
            print >>fout, "#PBS -l nodes=1:ppn=3"
            print >>fout, "#PBS -N %s" % pn
            print >>fout, "cd %s" % path1
            #print >>fout, 'debug/build_dist config.dk %(pi)d chr0'%locals() 
            #print >>fout, 'debug/alu_insert 1 config.dk %(pi)d '%locals()
            print >>fout, 'debug/alu_delete write_tmp1 config.dk %(pi)d'%locals() 
            print >>fout, 'debug/alu_delete write_tmp2 config.dk %(pi)d'%locals() 
            #print >>fout, 'debug/insert_pos 1 config.dk %(pi)d'%locals()
            print >>fout, 'debug/ins_del write_tmp1 config.dk %(pi)d'%locals() 
        pi += 1

if __name__ == '__main__':
    path1 = os.getcwd() + '/'
    pn_all = map(lambda x:x.strip(), file(path1 + 'inputs/PN_all').readlines())
    #print1(pn_all, path1+"q_dist/", path1)
    #print1(pn_all, path1+"q_ai/", path1)
    print1(pn_all, path1+"q_ad/", path1)
