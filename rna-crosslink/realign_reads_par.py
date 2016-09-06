import sys
import pysam
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from Bio import pairwise2
from multiprocessing import Pool

reffile = '/data/BOWTIE2_INDEXES/hg19.fa'
bamfile = 'work/head10k.hg19.bowtie.fix.bam'
seedsize = 40

if len(sys.argv) == 2:
    bamfile = sys.argv[1]
outfile = bamfile+'.realign'

def rev_com(s1):
    a = Seq(s1, generic_dna)
    return str(a.reverse_complement())

def align(s1, s2, step=9):
    ## fast align by a simple substring search 
    for i in xrange(step):
        a = i*(len(s1)/step)
        b = len(s1)-a-seedsize
        if b >= 0:
            r1 = s1[a:-b]
            if r1.find('-') < 0 and r1.find('[') < 0 and r1.find(']') < 0:
                A = s2.find(r1)
                S1 = s1
                S2 = s2
                if A >= 0: ## find a perfect match
                    if A > a:
                        S1 = '-'*(A-a)+S1
                    elif a > A:
                        S2 = '-'*(a-A)+S2
                    B = len(s2)-A-seedsize
                    if B > b:
                        S1 = S1+'-'*(B-b)
                    elif b > B:
                        S2 = S2+'-'*(b-B)
                    return [S1, S2]
    ## http://biopython.org/DIST/docs/api/Bio.pairwise2-module.html
    return pairwise2.align.localms(s1, s2, 1, -10, -10, -10)[0]

def comb_align(a1, a2, I=0, J=0):
    mul = align(a1[I].replace('-','['), a2[J].replace('-',']'))
    new = []
    for s in a1:
        n = ''
        j = 0
        for i in xrange(len(mul[0])):
            if mul[0][i] == '-':
                n += '-'
            else:
                n += s[j]
                j += 1
        new.append(n)
    for s in a2:
        n = ''
        j = 0
        for i in xrange(len(mul[1])):
            if mul[1][i] == '-':
                n += '-'
            else:
                n += s[j]
                j += 1
        new.append(n)
    return new

def mul_align(s1, g1, s2, g2, SPAN=-1):
    if SPAN < 0:
        a1 = align(s1, g1)[:2]
    else: ## is mapped
        a1 = ['-'*SPAN+s1+'-'*SPAN, g1]
    a2 = align(s1, g2)[:2]
    a3 = align(s2, g1)[:2]
    if SPAN < 0:
        a4 = align(s2, g2)[:2]
    else: ## is mapped
        a4 = ['-'*SPAN+s2+'-'*SPAN, g2]
    c1 = comb_align(a1, a2)[1:] ## [s1], g1, s1, g2
    c2 = comb_align(a3, a4)[1:] ## [s2], g1, s2, g2
    c3 = comb_align(c1, c2) ## g1, s1, g2 + g1, s2, g2
    return [c3[1], c3[0], c3[4], c3[5]]

def cmp_mul(s1, g1, g2):
    cc1 = 0
    cc2 = 0
    for i in xrange(len(s1)):
        if s1[i] != '-':
            if s1[i] == g1[i]: ## first align
                cc1 += 1
            elif s1[i] == g2[i]: ## second align
                cc2 += 1
    return cc1, cc2

ref = pysam.FastaFile(reffile)
sam = pysam.AlignmentFile(bamfile, 'rb')

data = []
chr_ref = list(sam.references)
last = sam.next()
for read in sam:
    if read.query_name != last.query_name: ## not a pair
        last = read
        continue
    f1 = bin(last.flag)
    f2 = bin(read.flag)
    if f1[-3] == '1' or f2[-3] == '1': ## not both mapped
        last = read
        continue
    assert f1[-7] == '1' and f2[-8] == '1'
    cig1 = last.cigar
    cig2 = read.cigar
    m1 = sum([j for i,j in cig1 if i==0])
    m2 = sum([j for i,j in cig2 if i==0])
    if m1 < seedsize or m2 < seedsize: ## bad mapping quality
        last = read
        continue
    r1 = chr_ref[last.reference_id]
    r2 = chr_ref[read.reference_id]
    p1 = last.reference_start
    p2 = read.reference_start
#    if r1 != r2 or abs(p1-p2) < 400:
    if r1 == r2:
        last = read
        continue
    s1 = last.query_sequence
    s2 = read.query_sequence
    if len(s1)-m1 < seedsize and len(s2)-m2 < seedsize: ## not enough overhang
        last = read
        continue
    if cig1[0][0] == 4:
        p1_clip = cig1[0][1]
    else:
        p1_clip = 0
    if cig2[0][0] == 4:
        p2_clip = cig2[0][1]
    else:
        p2_clip = 0
    p1A = p1-p1_clip
    p1B = p1-p1_clip+len(s1)
    p2A = p2-p2_clip
    p2B = p2-p2_clip+len(s2)
    data.append((read.query_name, r1, p1A, p1B, r2, p2A, p2B, f1, f2, s1, s2))
sam.close()

def run(par):
    SPAN = 50
    try:
        name, r1, p1A, p1B, r2, p2A, p2B, f1, f2, s1, s2 = par
        g1 = ref.fetch(r1, p1A-SPAN, p1B+SPAN)
        g2 = ref.fetch(r2, p2A-SPAN, p2B+SPAN)
        strand1 = '+'
        strand2 = '+'
        if not len(f1) >= 5 and f1[-5] == '1':
            strand1 = '-'
            s1 = rev_com(s1)
            g1 = rev_com(g1)
        if (len(f2) >= 5 and f2[-5] == '1'): ## RNA is the same as the 2nd read
            strand2 = '-'
            s2 = rev_com(s2)
            g2 = rev_com(g2)
        mu = mul_align(s1, g1, s2, g2, SPAN) ## multiple sequence alignment
        c1 = cmp_mul(mu[0], mu[1], mu[3])
        c2 = cmp_mul(mu[2], mu[3], mu[1])
        if (c1[0]>=seedsize and c1[1]>=seedsize and sum(c1)>=100) or (c2[0]>=seedsize and c2[1]>=seedsize and sum(c2)>=100):
            return strand1, strand2, c1, c2, mu
    except Exception,e:
        print e
        return

print 'There are', len(data), 'read pairs'

import time

pool = Pool()
time1 = time.time()
runlen = len(data) # min(len(data), 50000)
print 'We realign', runlen, 'for this run'
results = pool.map(run, data[:runlen])
diff = time.time() - time1
print 'We use', diff, 'seconds in this run'
print 'Run all samples will need', float(len(data))*diff/runlen/60/60, 'hours'
pool.close()
pool.join()

#results = [run(p) for p in data[:1]]

ref.close()

out = open(outfile, 'w')
for par, re in zip(data, results):
    if re is not None:
        name, r1, p1A, p1B, r2, p2A, p2B, f1, f2, s1, s2 = par
        strand1, strand2, c1, c2, mu = re
        out.write('>%s\n'%name)
        out.write('%s\n'%'\t'.join(['%s:%s-%s'%(r1, p1A+1, p1B), strand1, '%s:%s-%s'%(r2, p2A+1, p2B), strand2, str(c1), str(c2)]))
        out.write('%s\n%s\n%s\n%s\n'%tuple(mu))
        out.write('\n')
out.close()
