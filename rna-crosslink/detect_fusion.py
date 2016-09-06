import sys
import numpy as np

alignfile = 'work/head10k.hg19.bowtie.fix.bam.realign'
seedsize = 30

if len(sys.argv) == 2:
    alignfile = sys.argv[1]
outputfile = alignfile+'.fusion'

def match(r, g1, g2):
    m1 = np.zeros(len(r))
    m2 = np.zeros(len(r))
    for i in xrange(len(r)):
        if r[i] == '-':
            continue
        if r[i] == g1[i]:
            m1[i] = 1
        if r[i] == g2[i]:
            m2[i] = 1
    return m1, m2

def find(m1, m2):
    l = len(m1)
    c1 = np.cumsum(m1)
    c2 = np.cumsum(m2)
    if c1[l/2] > c2[l/2]:
## Case 1:
## m1 ||||| 
## m2      |||||
## c1 1234555555
## c2 0000012345
## df 1234543210
        df = c1-c2+m1+m2
        if df.sum() == 0:
            return 0
        return np.argmax(df)+1
    else:
## Case 2:
## m1      |||||
## m2 ||||| 
## c1 0000012345
## c2 1234555555
## df 1234543210
        df = c2-c1+m1+m2
        if df.sum() == 0:
            return 0
        return -(np.argmax(df)+1)

def location(f, g1, g2, loci1, loci2, s1, s2):
    ch1, p1ab = loci1.split(':')
    p1a, p1b = p1ab.split('-')
    ch2, p2ab = loci2.split(':')
    p2a, p2b = p2ab.split('-')
    assert f > 0
    k1 = sum([1 for i in g1[f:] if i != '-'])
    k2 = sum([1 for i in g2[:f] if i != '-'])
    out = []
    out.append(ch1)
    if s1 == '+':
        out.append(str(int(p1b)-k1))
    else:
        out.append(str(int(p1a)+k1))
    out.append(s1)
    out.append(ch2)
    if s2 == '+':
        out.append(str(int(p2a)+k2))
    else:
        out.append(str(int(p2b)-k2))
    out.append(s2)
    return out

def add_count(c, f, s, SPAN=30):
    S = s[(f-SPAN):(f+SPAN)] 
    for a in 'ACGT-':
        if a in c:
            count = c[a]
        else:
            count = [0 for i in xrange(len(S))]
        for i in xrange(len(S)):
            if S[i] == a:
                count[i] += 1
        c[a] = count

align = open(alignfile, 'r')
output = open(outputfile, 'w')
count = {}
cc = 0
while True:
    line = align.readline()
    if line == '':
        break
    if not line.startswith('>'):
        continue
    meta = align.readline().strip()
    loci1, strand1, loci2, strand2 = meta.split()[:4]
    r1 = align.readline().strip()
    g1 = align.readline().strip()
    r2 = align.readline().strip()
    g2 = align.readline().strip()
    m1a, m1b = match(r1, g1, g2)
    if ((1-m1a)*m1b).sum()>=seedsize and ((1-m1b)*m1a).sum()>=seedsize:
        f1 = find(m1a, m1b)
        if f1 > 0:
            fus = location(f1, g1, g2, loci1, loci2, strand1, strand2)
            output.write('\t'.join(fus+[line[1:-1], r1[:f1].lower()+r1[f1:]])+'\n')
            add_count(count, f1, r1)
        elif f1 < 0:
            fus = location(-f1, g2, g1, loci2, loci1, strand2, strand1)
            output.write('\t'.join(fus+[line[1:-1], r1[:-f1]+r1[-f1:].lower()])+'\n')
            add_count(count, -f1, r1)
    m2a, m2b = match(r2, g2, g1)
    if ((1-m2a)*m2b).sum()>=seedsize and ((1-m2b)*m2a).sum()>=seedsize:
        f2 = find(m2a, m2b)
        if f2 > 0:
            fus = location(f2, g2, g1, loci2, loci1, strand2, strand1)
            output.write('\t'.join(fus+[line[1:-1], r2[:f2]+r2[f2:].lower()])+'\n')
            add_count(count, f2, r2)
        elif f2 < 0:
            fus = location(-f2, g1, g2, loci1, loci2, strand1, strand2)
            output.write('\t'.join(fus+[line[1:-1], r2[:-f2].lower()+r2[-f2:]])+'\n')
            add_count(count, -f2, r2)
        cc += 1
align.close()
output.close()
print 'We have', cc, 'alignments'

import matplotlib.pyplot as plt
sumup = [i for i in count['A']]
for i in xrange(len(sumup)):
    for a in 'TCG':
        sumup[i] += count[a][i]
for a in 'ATCG':
    l = len(count[a])
    plt.plot(xrange(-l/2+1, l/2+1), [i/float(j) for i,j in zip(count[a], sumup)], label=a)
plt.xlabel('Fusion point')
plt.ylabel('Ratio')
plt.legend()
plt.savefig(alignfile+'.png')

