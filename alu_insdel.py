## build consensus running.
import sys, os

def read_config(fn, key):
    with open(fn) as fin:
        for line in fin:
            if line.startswith("#"):
                continue
            try:
                k, v = line.strip().split()[:2]
            except:
                continue
            if k == key:
                return v
    print 'ERROR', key, 'not exists at', fn
    sys.exit(0)

if __name__ == "__main__":
    if ( len(sys.argv) < 4):
        print 'usage: python %s config.dk chr1 q_cons'%(sys.argv[0])
        sys.exit(0)

    config_file = sys.argv[1]
    chrn = sys.argv[2]
    cmd_output_dir = sys.argv[3]
    
    if not cmd_output_dir.endswith('/'):
        cmd_output_dir += '/'

    cmd_output_fn = "tmpcons"   
    cmd_output_file_name = cmd_output_dir + cmd_output_fn 

    if os.path.exists(cmd_output_dir):
        print "warning, %s exists, might overwrite existing files"%cmd_output_dir
    else:
        os.mkdir(cmd_output_dir)

    cmd_file = read_config( config_file, "file_ins_cons") + chrn + ".pos"
    print 'rread', cmd_file
    print "qsub cmds will be output to %s*"%cmd_output_file_name


    s_bin = read_config( config_file, "bin_path") + 'debug/' + sys.argv[0].replace('.py','')
    cwd = os.getcwd()
    if cwd in config_file:
        s_config = config_file
    else:
        s_config = cwd + '/' + config_file
    s_pos = read_config( config_file, "file_ins_cons") + chrn + '_pos/'

    fout_sh = file(cmd_output_file_name + '.sh', 'w')
    i = 0
    with open(cmd_file) as fin:
        for line in fin:
            pos = line.strip()
            with open( '%s_%d'%(cmd_output_file_name, i), 'w') as fout:
                print >>fout, "#!/bin/sh "
                print >>fout, "#PBS -q normal"
                print >>fout, "#PBS -l nodes=1:ppn=1"
                print >>fout, "#PBS -l walltime=0:59:0"
                print >>fout, "#PBS -N %s_%d" % (cmd_output_fn, i)
                ###print >>fout, '/home/qianyuxx/faststorage/Alu/debug/alu_insdel /home/qianyuxx/faststorage/Alu/config.dk cons_reads_build '
                print >>fout, '%s %s cons_reads_build %s%s'%(s_bin, s_config, s_pos, pos)
            print >> fout_sh, 'qsub %s_%d'%(cmd_output_fn, i)
            i += 1    
    fout_sh.close()
