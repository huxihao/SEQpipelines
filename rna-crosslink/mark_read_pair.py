import sys, pysam

bam = 'work/head10k.hg19.bowtie.bam'

if len(sys.argv) == 2:
    bam = sys.argv[1]

inf = pysam.AlignmentFile(bam, 'rb')
ouf = pysam.AlignmentFile(bam[:-4]+'.fix.bam', 'wb', header=inf.header)
read1 = True
rname = ''
for a in inf:
    if rname == '':
        rname = a.query_name
    elif rname == a.query_name: ## the mete read
        read1 = False 
    if read1:
        a.flag += 64
    else:
        a.flag += 128
    ouf.write(a)
inf.close()
ouf.close()
