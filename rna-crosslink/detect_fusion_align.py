import sys
import numpy as np

alignfile = 'work/Control1.hg19.bowtie.fix.bam.realign'
seedsize = 30

if len(sys.argv) == 2:
    alignfile = sys.argv[1]
outputfile = alignfile+'.fusion'

def match(r, g1, g2):
    '''
        Find the number of matches to both genes
    '''
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
    '''
        Find out the fusion location
    '''
    l = len(m1)
    c1 = np.cumsum(m1)
    C1 = sum(m1)-c1
    c2 = np.cumsum(m2)
    C2 = sum(m2)-c2
    df1 = c1 + C2
## Case 1:
## m1 || || 
## m2      ||||
## c1 122344444
## c2 000001234
## C2 444443210
## df 566787654
    df2 = C1 + c2
## Case 2:
## m1      ||||
## m2 || || 
## c1 000001234
## C1 444443210
## c2 122344444
## df 566787654
    if df1.max() > df2.max():
        return 'Case1', [i+1 for i in xrange(l) if df1[i] == df1.max() and i>=5 and i<l-5]
    elif df2.max() > df1.max():
        return 'Case2', [i+1 for i in xrange(l) if df2[i] == df2.max() and i>=5 and i<l-5]
    else:
        return 'None', []

def location(f, g1, g2, loci1, loci2, s1, s2):
    '''
        Map the fusion location to genome
    '''
    ch1, p1ab = loci1.split(':')
    p1a, p1b = p1ab.split('-')
    ch2, p2ab = loci2.split(':')
    p2a, p2b = p2ab.split('-')
    k1 = sum([1 for i in g1[min(f):] if i != '-'])
    k2 = sum([1 for i in g2[:max(f)] if i != '-'])
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
count1 = {} ## reads with additional nucleotide
count2 = {} ## all fusion points

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
        c1, f1 = find(m1a, m1b)
        if len(f1) == 0:
            pass
        else:
            min_f1 = min(f1)
            max_f1 = max(f1)
            if c1 == 'Case1':
                G1, G2 = g1, g2
                fus = location(f1, g1, g2, loci1, loci2, strand1, strand2)
            elif c1 == 'Case2':
                G1, G2 = g2, g1
                fus = location(f1, g2, g1, loci2, loci1, strand2, strand1)
            output.write('\t'.join(fus+[line[1:-1]])+'\n')
            output.write(G1[:min_f1].lower()+G1[min_f1:]+'\n')
            output.write(r1[:min_f1].lower()+r1[min_f1:max_f1]+r1[max_f1:].lower()+'\n')
            findC1 = r1[(min_f1-1):(max_f1+1)].find('C')
            output.write(G2[:max_f1]+G2[max_f1:].lower()+'\n%s\n'%findC1)
            if findC1 < 0:
                findC1 = (min_f1+max_f1)/2
            else:
                findC1 += min_f1-1
            add_count(count2, findC1, r1)
            if max_f1 - min_f1 == 1 and r1[min_f1] != g1[min_f1] and r1[min_f1] != g1[min_f1]:
                add_count(count1, min_f1, r1)
                cc += 1
                continue ## no need to check read2
    m2a, m2b = match(r2, g2, g1)
    if ((1-m2a)*m2b).sum()>=seedsize and ((1-m2b)*m2a).sum()>=seedsize:
        c2, f2 = find(m2a, m2b)
        if len(f2) == 0:
            pass
        else:
            min_f2 = min(f2)
            max_f2 = max(f2)
            if c2 == 'Case1':
                G1, G2 = g2, g1
                fus = location(f2, g2, g1, loci2, loci1, strand2, strand1)
            elif c2 == 'Case2':
                G1, G2 = g1, g2
                fus = location(f2, g1, g2, loci1, loci2, strand1, strand2)
            output.write('\t'.join(fus+[line[1:-1]])+'\n')
            output.write(G1[:min_f2].lower()+G1[min_f2:]+'\n')
            output.write(r2[:min_f2].lower()+r2[min_f2:max_f2]+r2[max_f2:].lower()+'\n')
            findC2 = r2[(min_f2-1):(max_f2+1)].find('C')
            output.write(G2[:max_f2]+G2[max_f2:].lower()+'\n%s\n'%findC2)
            if findC2 < 0:
                findC2 = (min_f2+max_f2)/2
            else:
                findC2 += min_f2-1
            add_count(count2, findC2, r2)
            if max_f2 - min_f2 == 1 and r2[min_f2] != g1[min_f2] and r2[min_f2] != g2[min_f2]:
                add_count(count1, min_f2, r2)
                cc += 1
align.close()
output.close()
print 'We have', cc, 'alignments'

def plot_ratio(count, fname):
    import matplotlib.pyplot as plt
    sumup = [i for i in count['A']]
    for i in xrange(len(sumup)):
        for a in 'TCG':
            sumup[i] += count[a][i]
    for a in 'ATCG':
        l = len(count[a])
        plt.plot(xrange(-l/2, l/2), [i/float(j) for i,j in zip(count[a], sumup)], label=a)
    plt.title('Plot for %s reads'%max(sumup))
    plt.ylim([0, 0.8])
    plt.xlabel('Fusion point')
    plt.ylabel('Ratio')
    plt.legend()
    plt.savefig(fname+'.png')
    plt.clf()

plot_ratio(count1, alignfile+'.addition')
plot_ratio(count2, alignfile+'.findC')

