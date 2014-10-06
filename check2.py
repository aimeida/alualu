## new output format, vcf only
import sys
from scipy.stats import itemfreq
import itertools

def family_member(gid):
    if gid == "01":
        return 'F'
    elif gid == "02":
        return 'M'
    return 'C'

class VCF_PARSER():
    def _parse_trio_group(self):
        for pn in self.pns:
            gname, gid = pn.split('-')
            if gname not in self.trio_group:
                self.trio_group[gname] = {}
            self.trio_group[gname][family_member(gid)] = pn
            self.seqs[pn] = ''

    def __init__(self, fn_vcf, fn_pn_used):
        with open(fn_vcf) as fin:
            for line in fin:
                if line.startswith('#CHROM'):
                    header = line.strip().split('\t')
                elif not line.startswith('#'):
                    break
        self.fn_vcf = fn_vcf
        self.pns = map(lambda x:x.strip(), file(fn_pn_used).readlines())
        self.idx = map(lambda p: header.index(p) if p in header else -1, self.pns)
        self.header = header
        self.trio_group = {}
        self.chr_pos = []
        self.seqs = {}
        self._parse_trio_group()
        

    def read_content(self, qual_th):
        id1 = self.header.index('FILTER')
        with open(self.fn_vcf) as fin:
            for line in fin:
                if line.startswith('#'):
                    continue
                tmp = line.strip().split('\t')
                if not tmp[id1] in qual_th:
                    continue
                self.chr_pos.append( tmp[0]+'_'+tmp[1] )
                for i, j in zip(self.pns, self.idx):
                    if j < 0:
                        continue
                    self.seqs[i] += tmp[j].split(':')[0]
    
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
        if '.' in [pn1[i], pn2[i], pn3[i]]:
            continue
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
##        if pos[i] == "chr1_11096562":
##            print from1, from2, pn1[i], pn2[i], pn3[i]

    conflict_rate = len(conflict_pos) / float(npos)
    print '%d positions OK, confliction rate: %.1f %% for %d positions' %(npos - len(conflict_pos), conflict_rate * 100, npos)
    if verbose:
        print '\n'.join(conflict_pos)


def read_vcf(fn, qual_th, fn_pn_used):
    vcf_fh = VCF_PARSER(fn_vcf, fn_pn_used)
    vcf_fh.read_content(qual_th)
    #print len(vcf_fh.seqs), len(vcf_fh.chr_pos)
    #print vcf_fh.seqs.values()[0]
    for gn, v1 in vcf_fh.trio_group.items():
        if verbose:
            #if gn != '1006':
            if gn != '1473':
                continue                                                                                                             
        if len(v1) != 3:
            continue
        pn_father = v1['F']
        pn_mother = v1['M']
        pn_child = v1['C']
        print 'check family ', gn
        one_family(vcf_fh.seqs[pn_father], vcf_fh.seqs[pn_mother], vcf_fh.seqs[pn_child], vcf_fh.chr_pos, allow_denovo, verbose)


def hwe_txt(fn_vcf, fn_pn_used, fn_txt, qual_th):
    vcf_fh = VCF_PARSER(fn_vcf, fn_pn_used)
    id1 = vcf_fh.header.index('FILTER')
    id2 = vcf_fh.header.index('FORMAT') + 1
    #id3 = vcf_fh.header.index('INFO')
    fout = file(fn_txt,'w')
    print >>fout, '#CHROM POS', ' '.join(vcf_fh.pns)
    with open(fn_vcf) as fin:
        for line in fin:
            if line.startswith('#'):
                continue
            tmp = line.strip().split('\t')
            ###af = tmp[id3].split('AF=')[1]
            if tmp[id1] in qual_th:
                genos = map(lambda x:x.split(':')[0], tmp[id2:])
                if '.' not in genos:
                    print >>fout, tmp[0], tmp[1], ' '.join(genos)

if __name__ == "__main__":
    
    opt = sys.argv[1]

    if opt == 'i':
        print "next step:  sh check.sh 0 1006 chr1 44059297"
        file_pn_used = '/home/qianyuxx/faststorage/AluDK/outputs/insert_alu1/pn_used'
        fn_vcf = '/home/qianyuxx/faststorage/AluDK/outputs/insert_alu1/fixed_delete0/29.vcf'
                
    elif opt == 'd':
        print "next step:  sh check.sh 1 1006 chr1 44059297"
        file_pn_used = '/home/qianyuxx/faststorage/AluDK/inputs/PN_all'
        fn_vcf = '/home/qianyuxx/faststorage/AluDK/outputs/delete_alu0/30.vcf' 

    allow_denovo = True 
    #allow_denovo = False
    #verbose = False
    verbose = True
    qual_th = ['PASS']
    #qual_th = ['PASS', 'BreakpointOneside']

    print 'checking vcf file of', fn_vcf
    print 'allow denovo', allow_denovo
    read_vcf(fn_vcf, qual_th, file_pn_used)
    hwe_txt(fn_vcf, file_pn_used, fn_vcf.replace('.vcf', '.txt'), qual_th) 
