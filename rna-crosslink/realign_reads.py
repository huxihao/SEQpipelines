import sys
import pysam
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from Bio import pairwise2

reffile = '/data/BOWTIE2_INDEXES/hg19.fa'
bamfile = 'work/head10k.hg19.bowtie.fix.bam'
seedsize = 35

if len(sys.argv) == 2:
    bamfile = sys.argv[1]
outfile = bamfile+'.realign'

def rev_com(s1):
    a = Seq(s1, generic_dna)
    return str(a.reverse_complement())

def align(s1, s2):
    align = pairwise2.align.localms(s1, s2, 1, -100, -10, -100)
#    print pairwise2.format_alignment(*align[0])
    return align[0]

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

def mul_align(s1, g1, s2, g2):
#    a1 = align(s1, g1)[:2]
    a1 = [s1, g1]
    a2 = align(s1, g2)[:2]
    a3 = align(s2, g1)[:2]
#    a4 = align(s2, g2)[:2]
    a4 = [s2, g2]
    c1 = comb_align(a1, a2)[1:] ## [s1], g1, s1, g2
    c2 = comb_align(a3, a4)[1:] ## [s2], g1, s2, g2
    c3 = comb_align(c1, c2) ## g1, s1, g2 + g1, s2, g2
    return [c3[1], c3[0], c3[4], c3[5]]

def mul_align_old(s1, g1, g2):
    a1 = align(s1, g1)
    a2 = align(s1, g2)
    mul = align(a1[0].replace('-','['), a2[0].replace('-',']'))
    mul_s1 = mul[0].replace('[','-')
    mul_g1 = ''
    mul_g2 = ''
    j = 0
    k = 0
    for i in xrange(len(mul[0])):
        m1 = mul[0][i]
        m2 = mul[1][i]
        if m1 == '-':
            mul_g1 += '-'
        else:
            if m1 == a1[1][j]:
                mul_g1 += m1
            else:
                mul_g1 += a1[1][j].lower()
            j += 1
        if m2 == '-':
            mul_g2 += '-'
        else:
            if m2 == a2[1][k]:
                mul_g2 += m2
            else:
                mul_g2 += a2[1][k].lower()
            k += 1
    return mul_s1, mul_g1, mul_g2

def cmp_mul(s1, g1, g2):
    cc1 = 0
    cc2 = 0
    for i in xrange(len(s1)):
        if s1[i] != '-':
            if s1[i] == g1[i]:
                cc1 += 1
            if s1[i] == g1[i] or s1[i] == g2[i]:
                cc2 += 1
    return cc1, cc2

ref = pysam.FastaFile(reffile)
sam = pysam.AlignmentFile(bamfile, 'rb')
out = open(outfile, 'w')

chr_ref = list(sam.references)
cc1 = 0
cc2 = 0
last = sam.next()
for read in sam:
    if read.query_name != last.query_name: ## not a pair
        last = read
        continue
    cc1 += 1
    f1 = bin(last.flag)
    f2 = bin(read.flag)
    if f1[-3] == '1' or f2[-3] == '1': ## not both mapped
        last = read
        continue
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
    if r1 == r2: ## skip intra
        last = read
        continue
    s1 = last.query_sequence
    s2 = read.query_sequence
    if cig1[0][0] == 4:
        p1_clip = cig1[0][1]
    else:
        p1_clip = 0
    if cig1[0][0] == 4:
        p2_clip = cig2[0][1]
    else:
        p2_clip = 0
    p1A = p1-p1_clip
    p1B = p1-p1_clip+len(s1)
    p2A = p2-p2_clip
    p2B = p2-p2_clip+len(s2)
    try:
        g1 = ref.fetch(r1, p1A, p1B)
        g2 = ref.fetch(r2, p2A, p2B)
    except: ## fail to fetch the genome sequence
        last = read
        continue
    strand1 = '+'
    strand2 = '+'
    if len(f1) >= 5 and f1[-5] == '1':
        strand1 = '-'
        s1 = rev_com(s1)
        g1 = rev_com(g1)
    if len(f2) >= 5 and f2[-5] == '1':
        pass
    else: ## read 2 is reversed
        strand2 = '-'
        s2 = rev_com(s2)
        g2 = rev_com(g2)
    mu = mul_align(s1, g1, s2, g2) ## multiple sequence alignment
    c1 = cmp_mul(mu[0], mu[1], mu[3])
    c2 = cmp_mul(mu[2], mu[3], mu[1])
    if c1[1] >= len(s1)-10 and c2[1] >= len(s2)-10 and (c1[1]-c1[0] >= seedsize or c2[1]-c2[0] >= seedsize):
        cc2 += 1
#        print '[%s, %s]'%(cc1, cc2), r1, p1, strand1, r2, p2, strand2, m1, c1, m2, c2
        out.write('>%s\n'%read.query_name)
        out.write('%s\n'%'\t'.join(['%s:%s-%s'%(r1, p1A+1, p1B), strand1, '%s:%s-%s'%(r2, p2A+1, p2B), strand2, str(c1), str(c2)]))
        out.write('%s\n%s\n%s\n%s\n'%tuple(mu))
        out.write('\n')
    last = read
ref.close()
sam.close()
out.close()
