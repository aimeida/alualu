import sys, os

fn_all = '/nfs_mount/bioinfo/users/yuq/work/Alu/inputs/PN_all'
pn_all = map(lambda x:x.strip(), file(fn_all).readlines())
idx_name = dict( [ (i,j) for i,j in enumerate(file(fn_all).read().split())] )
#pn_idx = map(int, file('/nfs_mount/bioinfo/users/yuq/work/Alu/inputs/PN_20').readlines())
pn_idx = range(len(pn_all))
#pn_idx = range(1000)

#fn_qsub = 'tmp1.txt'
#fn_qsub = 'qsub1.txt'
fn_qsub = 'qsub0.txt'

path1 = "/nfs_mount/bioinfo/users/yuq/work/Alu/stat/alu/"
path2 = "/nfs_mount/bioinfo/users/yuq/work/Alu/outputs/insert_len/count/"
path3 = "/nfs_mount/bioinfo/users/yuq/work/Alu/outputs/insert_len/prob/"
#path0 = "/nfs_mount/bioinfo/users/yuq/work/Alu/stat/alu/rp"  ## for deletion detection
path0 = "/nfs_mount/bioinfo/users/yuq/work/Alu/stat/alu/ip"   ## for insertion detection

# print '%(path1)sdebug/alu_insert 1 %(path1)sconfig.properties 0'%locals()
# sys.exit() 

pi = 0
for i in pn_idx:
    if pi % 500 == 0:
        if pi > 0:
            fout.close()
        path0_ext = path0+'%d/'%(i/500)
        if not os.path.exists(path0_ext):
            os.mkdir(path0_ext)
        fout = file(path0_ext + fn_qsub, 'w')
    pi += 1
    #print >>fout, '%(path1)sdebug/build_dist %(path1)sconfig.properties %(i)d chr0'%locals()
    #print >>fout, '%(path1)sdebug/alu_delete 1 %(path1)sconfig.properties %(i)d chr0'%locals()
    #print >>fout, '%(path1)sdebug/alu_insert 1 %(path1)sconfig.properties %(i)d'%locals()
    print >>fout, '%(path1)sdebug/insert_pos 1 %(path1)sconfig.properties %(i)d'%locals()
