import sys
import math
import itertools

def read_geno(v):
    vs = map(lambda x: math.exp(-float(x)/10), v.split(','))
    vs = map(lambda x: '%.2f'%(round(x*100)/100), vs)
    return vs

def call_geno(v):
    v = map(float, v.split(','))
    return str(v.index(min(v)))

def get_idx(fn, pns):
    with open(fn) as fin:
        for line in fin:
            if not line.startswith('#CHROM'):
                continue
            header = line.strip().split()
            idx = map(lambda p: header.index(p) if p in header else -1, pns)
            return idx

def parse_lines(f_in, fout, check_idx):
    with open(f_in) as fin:
        for line in fin:
            if line.startswith('#'):
                continue
            tmp = line.strip().split()
            vs =  map(call_geno, [tmp[i] for i in check_idx]) 
            print >>fout, tmp[0], tmp[1], ' '.join(vs)                
            
def vcf_to_text(pns, idx, f_in, f_out):
    idx_pn = dict([(j,i) for i,j in zip(pns, idx) if j >= 0])
    check_idx = sorted(idx_pn.keys())
    with open(f_out,'w') as fout:
        print >>fout, '#CHROM  POS', ' '.join([idx_pn[k] for k in check_idx])
        parse_lines(f_in, fout, check_idx)

def read_pn(f_vcf):
    with open(f_vcf) as fin:
        for line in fin:
            if line.startswith('#CHR'):
                tmp = line.strip().split()
                coln = tmp.index('FORMAT') + 1
                return tmp[coln:]

def family_member(gid):
    if gid == "01":
        return 'F'
    elif gid == "02":
        return 'M'
    return 'C'
        
def parse_trio_group(f_vcf):
    pns = read_pn(f_vcf)
    trio_group = {}  
    for pn in pns:
        gname, gid = pn.split('-')
        if gname not in trio_group:
            trio_group[gname] = {}
        trio_group[gname][family_member(gid)] = pn
    return trio_group


def filter_chisq(f_llh, f_input, f_output, offset=1):
    fin2 = file(f_input)
    fout = file(f_output,'w')
    fout.write(fin2.readline())

    with open(f_llh) as fin:
        fin.readline()
        for line in fin:
            line2 = fin2.readline()
            if float(line.strip().split()[-1]) >= 1.92: ### qchisq(0.95, df=1)  [1] 3.841459
                fout.write(line2)
    fin2.close()
    fout.close()

def read_vcftxt(f_vcftxt):
    chr_pos = []
    seqs = {}
    with open(f_vcftxt) as fin:
        pns = fin.readline().strip().split()[2:]
        for pn in pns:
            seqs[pn] = ''
        for line in fin:
            tmp = line.strip().split()
            chr_pos.append( tmp[0]+'_'+tmp[1] )
            for i,j in zip(pns, tmp[2:]):
                seqs[i] += j
    return chr_pos, seqs

def one_family(pn1, pn2, pn3, pos, allow_denovo, verbose):
    def p2c(parent):
        if parent == '2':
            return [1]
        elif parent == '0':
            return [0]
        elif parent == '1':
            return [0,1]
        else:
            print "ERROR!"
    nlen = len(pn1)
    npos = 0
    conflict_pos = []
    for i in range(nlen):
        from1 = p2c(pn1[i])
        from2 = p2c(pn2[i])
        if (pn1[i] == '0' and pn2[i] == '0' and pn3[i] == '0'):
            continue

        npos += 1
        if int(pn3[i]) not in [sum(x) for x in itertools.product(from1, from2)]:
            if allow_denovo: 
                if not (pn1[i] == '0' and pn2[i] == '0' and pn3[i] == '1'):
                    conflict_pos.append(pos[i])
            else:
                conflict_pos.append(pos[i])

    conflict_rate = len(conflict_pos) / float(npos)
    print '%d positions OK, confliction rate: %.1f %% for %d positions' %(npos - len(conflict_pos), conflict_rate * 100, npos)
    if verbose:
        print '\n'.join(conflict_pos)


if __name__ == "__main__":
    
    opt = sys.argv[1]
    
    allow_denovo = True 
    #allow_denovo = False
    verbose = False

    if opt == '1':
        file_pn_used = '/home/qianyuxx/faststorage/AluDK/outputs/insert_alu1/pn_used'
        f_llh = '/home/qianyuxx/faststorage/AluDK/outputs/insert_alu1/fixed_delete0/29.pos'
        f_vcf = '/home/qianyuxx/faststorage/AluDK/outputs/insert_alu1/fixed_delete0/29.vcf'
        
    elif opt == '1b':
        file_pn_used = '/home/qianyuxx/faststorage/AluDK/outputs/insert_alu1/pn_used'
        f_llh = '/home/qianyuxx/faststorage/AluDK/outputs/insert_alu1/fixed_delete0_backup/29.pos'
        f_vcf = '/home/qianyuxx/faststorage/AluDK/outputs/insert_alu1/fixed_delete0_backup/29.vcf'
        
    elif opt == '2':
        file_pn_used = '/home/qianyuxx/faststorage/AluDK/inputs/PN_all'
        f_llh = '/home/qianyuxx/faststorage/AluDK/outputs/delete_alu0/30.pos' 
        f_vcf = '/home/qianyuxx/faststorage/AluDK/outputs/delete_alu0/30.vcf' 

    elif opt == '3':
        file_pn_used = '/home/qianyuxx/faststorage/AluDK/inputs/PN_all'
        f_llh = '/home/qianyuxx/faststorage/AluDK/outputs/_Maj_delete_alu0/30.pos' 
        f_vcf = '/home/qianyuxx/faststorage/AluDK/outputs/_Maj_delete_alu0/30.vcf' 

        
    print 'checking vcf file of', f_vcf
    print 'allow denovo', allow_denovo

    pn_used = map(lambda x:x.strip(), file(file_pn_used).readlines())
    f_vcftxt = f_vcf + '.txt'
    f_vcftxt2 = f_vcftxt + '.filter'

    idx = get_idx(f_vcf, pn_used)
    # vcf_to_text(pn_used, idx, f_vcf, f_vcftxt)        
    # filter_chisq(f_llh, f_vcftxt, f_vcftxt2) 
    
    npos = len(file(f_vcftxt2).readlines()) - 1
    print npos, 'positions considered'

    trio_group= parse_trio_group(f_vcf)
    chr_pos, seqs = read_vcftxt(f_vcftxt2)
    for gn, v1 in trio_group.items():
        
        verbose = True
        if gn != '1006': 
            continue
        
        if len(v1) != 3:
            continue
        pn_father = v1['F']
        pn_mother = v1['M']
        pn_child = v1['C']
        print 'check family ', gn
        one_family(seqs[pn_father], seqs[pn_mother], seqs[pn_child], chr_pos, allow_denovo, verbose)
        
