#! /usr/bin/python
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

if __name__ == "__main__":
    
    opt = sys.argv[1]
    
    if opt == 'decode':
        file_pn_used = '/nfs/gpfs/data/pbstmp/bjarnih/Alu/140622/PN_used_del'
        f_llh = '/nfs/gpfs/data/pbstmp/bjarnih/Alu/140622/delete_alu0/2649.pos'
        f_vcf = '/nfs/gpfs/data/pbstmp/bjarnih/Alu/140622/delete_alu0/2649.vcf'
        
    print 'checking vcf file of', f_vcf
    print 'allow denovo', allow_denovo

    pn_used = map(lambda x:x.strip(), file(file_pn_used).readlines())
    f_vcftxt = f_vcf + '.txt'
    f_vcftxt2 = f_vcftxt + '.filter'

    idx = get_idx(f_vcf, pn_used)
    vcf_to_text(pn_used, idx, f_vcf, f_vcftxt)        
    filter_chisq(f_llh, f_vcftxt, f_vcftxt2) 
    
