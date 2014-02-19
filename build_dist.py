import sys, os
import numpy as np

def read_pn():
    fn = '/nfs_mount/bioinfo/users/yuq/work/Alu/inputs/PN_noP'
    return map(lambda x:x.strip(), file(fn).readlines()) 

def update_count(fn, dic):
    with open(fn) as fin:
        for line in fin:
            a, b = map(int, line.split())
            try:
                dic[a] += b
            except:
                dic[a] = b

def parse_file(fn):
    counts = 0 
    nsum = 0
    with open(fn) as fin:
        for line in fin:
            a, b = map(int, line.split())
            counts += b
            nsum += (a*b)
    return nsum/float(counts), counts

if __name__ == "__main__":
    opt = sys.argv[1]
    if opt == "1":
        with open(sys.argv[0].replace('.py','.txt'), 'w') as fout:
            #for i in range(100):
            for i in range(20):
                print >>fout, "/nfs_mount/bioinfo/users/yuq/work/Alu/stat/alu/debug/build_dist /nfs_mount/bioinfo/users/yuq/work/Alu/inputs/PN_noP /nfs/gpfs/data/Results/GWS/ /nfs_mount/bioinfo/users/yuq/work/Alu/outputs/insert_len/ chr13", i
            
    elif opt == "2": ## list Read groups 
        rg_dict = {}
        path = '/nfs_mount/bioinfo/users/yuq/work/Alu/outputs/insert_len/'
        for fn in os.listdir(path):
            rg = fn.split('chr1.')[1]
            try:
                rg_dict[rg] += 1
            except:
                rg_dict[rg] = 1
        print len(rg_dict), sum(rg_dict.values()) ## read group is not share between diff PN
    
    elif opt == "3": ## counts of reading groups for each PN
        chrx = sys.argv[2] ##"chr1"
        path = '/nfs_mount/bioinfo/users/yuq/work/Alu/outputs/insert_len/'
        n_pn = 20
        f_mean = ''
        for pn in read_pn()[:n_pn]:
            v_mean, v_count = [], []
            for fn in os.listdir(path):
                if pn in fn and '.'+chrx+'.' in fn:
                    _v =  parse_file(path + fn)                   
                    v_mean.append(_v[0])
                    v_count.append(_v[1])
            if v_mean:
                if not f_mean:
                    f_mean = file('/nfs_mount/bioinfo/users/yuq/work/Alu/outputs/grocery/RG_mean.'+chrx, 'w')     
                    f_count =  file('/nfs_mount/bioinfo/users/yuq/work/Alu/outputs/grocery/RG_count.'+chrx, 'w')
                print >> f_mean, pn, ' '.join(['%.1f'%i for i in v_mean])
                print >> f_count, pn, sum(v_count), ': ',' '.join(['%d'%i for i in v_count])
        if f_mean:
            f_mean.close()
            f_count.close()
    
    elif opt == "4": ## combine rg for chr1
        chrx = sys.argv[2] ##"chr1"
        path = '/nfs_mount/bioinfo/users/yuq/work/Alu/outputs/insert_len/%s/'%chrx        
        n_pn = 20
        for pn in read_pn()[:n_pn]:
            fn_new = ''
            v_count = {}
            for fn in os.listdir(path):
                if pn in fn and '.'+chrx+'.' in fn:
                    rg = fn.split('chr1.')[1]
                    if not fn_new:
                        fn_new = path + fn.replace('.'+rg, '')
                    if rg:
                        update_count(path+fn, v_count)
            print 'output into', fn_new
            with open(fn_new, 'w') as fout:
                for key in sorted(v_count.keys()):
                    print >> fout, key, v_count[key]
    elif opt == "5": ## print cmd, not used any more
        path1 = '/nfs_mount/bioinfo/users/yuq/work/Alu/outputs/insert_len/chr1/'
        binpath = '/nfs_mount/bioinfo/users/yuq/work/Alu/stat/alu/debug/build_dist'
        rpath = '/nfs_mount/bioinfo/users/yuq/work/Alu/cmds/InsertLenDist.R'
        n_pn = 20
        for pn in read_pn()[:n_pn]:
            f_in = path1 + pn + '.count.chr1'
            print '%s %s %s.100_1000_5 100 1000 5'%(binpath, f_in, f_in) ## pdf_file
            print 'Rscript %s %s.100_1000_5'%(rpath, f_in)              ## plot log(p)

    elif opt == "6":  ## build complete pdf for each RG
        path1 = '/nfs_mount/bioinfo/users/yuq/work/Alu/outputs/insert_len/20_chr/'
        path2 = '/nfs_mount/bioinfo/users/yuq/work/Alu/outputs/insert_len/prob_20_chr/'
        binpath = '/nfs_mount/bioinfo/users/yuq/work/Alu/stat/alu/debug/build_dist'
        rpath = '/nfs_mount/bioinfo/users/yuq/work/Alu/cmds/InsertLenDist.R'
        fns = os.listdir(path1)
        for fn in fns:
            print '%s %s %s.100_1000_5 100 1000 5'%(binpath, path1 + fn, path2 + fn) ## pdf_file
            print 'Rscript %s %s.100_1000_5'%(rpath, path2 + fn)              ## plot log(p)
