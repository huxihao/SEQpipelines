import os
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--datapath', type=str, default='test/', help='Data path')
parser.add_argument('--workpath', type=str, default='work/', help='Work path')
parser.add_argument('--genome', type=str, default='hg19', help='Reference genome')
parser.add_argument('--threads', type=int, default=10, help='Number of threads')
parser.add_argument('--server', type=str, default='../server.txt', help='File for saving the server information')

args = parser.parse_args()
print args

genome = args.genome
thread = args.threads
datapath = args.datapath
workpath = args.workpath
server = {}
with open(args.server, 'r') as tmp:
    for line in tmp:
        ele = line.split()
        server[ele[0]] = ele[1]

os.environ['PATH'] = os.environ['PATH']+':'+os.path.abspath('bin')
os.environ['BOWTIE2_INDEXES'] = '/data/BOWTIE2_INDEXES/'

htmlmain = open(workpath+'index.html', 'w')

with open('../styles.css', 'r') as temp:
    lines = temp.readlines()
htmlmain.write('''
<head>
<style type="text/css">
%s
</style>
'''%''.join(lines))

dataset = []
for p in os.listdir(datapath):
    if os.path.isfile(datapath+p):
        continue
    f = [i for i in os.listdir(datapath+p) if i.endswith('.fastq.gz')]
    f.sort()
    print p, f
    dataset.append((p, datapath+p+'/'+f[0], datapath+p+'/'+f[1]))
dataset.sort()

htmlmain.write('''
<h1>Report for Cross-linking RNA-seq</h1>
<h2>Data Files</h2>
<table>
<tr><th>Name</th><th>Read 1</th><th>Read 2</th></tr>
''')
for d, f1, f2 in dataset:
    htmlmain.write('<tr><th>%s</th><td>%s</td><td>%s</td></tr>'%(d,f1,f2))

htmlmain.write('''
</table>        
        
<h2>Data Quality</h2>
<table>
<tr><th>Name</th><th>Fastx output</th><th>Fastqc output</th></tr>
 ''')
for d, f1, f2 in dataset:
    print '> Check', d, 'to', genome
    out = '%s.%s.fastx.txt'%(d, genome)
    if not os.path.exists(workpath+out):
        os.system('zcat '+f1+' | fastx_quality_stats > '+workpath+out)
        os.system('zcat '+f2+' | fastx_quality_stats >>'+workpath+out)
    f1name = (f1.split('/')[-1]).split('.')[0]+'_fastqc'
    if not os.path.exists(workpath+f1name):
        os.system('fastqc '+f1+' -o '+workpath)
    f2name = (f2.split('/')[-1]).split('.')[0]+'_fastqc'
    if not os.path.exists(workpath+f2name):
        os.system('fastqc '+f2+' -o '+workpath)
    htmlmain.write('<tr><th>%s</th><td><a href="%s">%s</a></td><td><a href="%s/fastqc_report.html" target="_blank">Read 1</a> <a href="%s/fastqc_report.html" target="_blank">Read 2</a></td></tr>'%(d, out, out, f1name, f2name))
htmlmain.write('''
</table>        

<h2>Map reads by Bowtie</h2>
<table>
<tr><th>Name</th><th>Command</th></tr>
''')
for d, f1, f2 in dataset:
    print '> Bowtie', d, 'to', genome
    out = workpath+'%s.%s.bowtie.bam'%(d, genome)
    cmd = 'zcat '+f1+' '+f2+' | '+'bowtie2 '+genome+' -p '+str(thread)+' -k 1 --local --reorder - | samtools view -ubS - -o '+out
    if not os.path.exists(out):
        os.system(cmd)
    htmlmain.write('<tr><th>%s</th><td>%s</td></tr>'%(d, cmd))
htmlmain.write('''
</table>

<h2>Fix pair-end reads</h2>
<table>
<tr><th>Name</th><th>Total</th><th>Mapped</th><th>Paired</th><th>Log</th></tr>
''')
for d, f1, f2 in dataset:
    print '> Fix', d
    pre = workpath+'%s.%s.bowtie'%(d, genome)
    if not os.path.exists(pre+'.fix.bam.log'):
        os.system('python mark_read_pair.py '+pre+'.bam')
        os.system('samtools sort -@ '+str(thread)+' -n '+pre+'.fix.bam > '+pre+'.srt.bam') ## sort by name
        os.system('samtools fixmate '+pre+'.srt.bam '+pre+'.fix.bam')
        os.system('samtools flagstat '+pre+'.fix.bam > '+pre+'.fix.bam.log')
        os.system('samtools sort -@ '+str(thread)+' '+pre+'.fix.bam | samtools view -buh -F 4 - | samtools view -bh -F 8 - > '+pre+'.srt.bam') ## sort and index
        os.system('samtools index '+pre+'.srt.bam')
    with open(pre+'.fix.bam.log', 'r') as tmp:
        lines = tmp.readlines()
    htmlmain.write('<tr><th>%s</th><td>%s</td><td>%s</td><td>%s</td><td><a href="%s.fix.bam.log">Open</a></td></tr>'%(d, lines[0].split()[0], lines[4].split()[0], lines[9].split()[0], pre[len(workpath):]))
htmlmain.write('''
</table>

<h2>Genome browser tracks</h2>
<p>Copy the following tracks:</p>
<table><tr><td>''')
for d, f1, f2 in dataset:
    print '> Output Genome Browser tracks'
    out = '%s.%s.bowtie.srt.bam'%(d, genome)
    htmlmain.write('track type=bam pairEndsByName=. db=%s name="%s" visibility=2 bamColorMode=gray bigDataUrl="%s/%s/rna-crosslink/%s"<br/>'%(genome, d, server['address'], server['home'], out))
htmlmain.write('''</td></tr></table>
<p>into <a href="http://genome.ucsc.edu/cgi-bin/hgCustom?db=%s" target="_black">the UCSC genome browser</a>, or use <a href="%s/internal/browser/cgi-bin/hgCustom?db=%s" target="_black">the local genome browser</a> (user name: xuelab; password: neibu)</p>

<h2>Read pair summary</h2>
<table>
<tr><th>Name</th><th>All pairs</th><th>Within 400bp</th><th>More than 400bp</th><th>Different chromosomes</th><th>Histogram</th></tr>
'''%(genome, server['address'], genome))
for d, f1, f2 in dataset:
    print '> Summary', d
    pre = workpath+'%s.%s.bowtie.fix'%(d, genome)
    cmd = 'bedtools bamtobed -i '+pre+r'.bam -bedpe -mate1 -ed | grep -v [.] | gzip -cf > '+pre+'.gz'
    if not os.path.exists(pre+'.gz'):
        os.system(cmd)
    import numpy as np
    if not os.path.exists(pre+'.npy'):
        data = []
        import gzip
        tmp = gzip.open(pre+'.gz', 'rb')
#        known = set()
        for line in tmp:
            ele = line.split()
            ch1, p1a, p1b, ch2, p2a, p2b, rname, mismatch, st1, st2 = ele
            if ch1 != '.' and ch2 != '.':
                if ch1 == ch2:
                    all_d = [int(p1a), int(p1b), int(p2a), int(p2b)]
                    dis = max(all_d) - min(all_d) ## fragment size
                else:
                    dis = 0
                if st1 == st2.strip():
                    dis = -dis
#                if (ch1, p1a, p1b, ch2, p2a, p2b) not in known:
                data.append(dis)
#                    known.add((ch1, p1a, p1b, ch2, p2a, p2b))
        tmp.close()
        np.save(pre+'.npy', np.array(data))
    import matplotlib.pyplot as plt
    dd = np.load(pre+'.npy')
    abdd = np.abs(dd)
    cc_total = np.sum(dd==dd)
    cc_less400 = np.sum(np.logical_and(0<abdd, abdd<=400))
    cc_more400 = np.sum(abdd>400)
    cc_diffchr = np.sum(abdd==0)
    case1a = abdd[(dd>0)&(0<abdd)&(abdd<1000)]
    case1b = abdd[(dd<0)&(0<abdd)&(abdd<1000)]
    case2a = abdd[(dd>0)&(1000<abdd)&(abdd<10000)]
    case2b = abdd[(dd<0)&(1000<abdd)&(abdd<10000)]
    case3a = abdd[(dd>0)&(10000<abdd)&(abdd<100000)]
    case3b = abdd[(dd<0)&(10000<abdd)&(abdd<100000)]
    plt.subplot(3,1,1)
    plt.hist([case1a, case1b] if len(case1b)>0 else case1a, 100, stacked=True)
    plt.subplot(3,1,2)
    plt.hist([case2a, case2b] if len(case2b)>0 else case2a, 100, stacked=True)
    plt.subplot(3,1,3)
    plt.hist([case3a, case3b] if len(case3b)>0 else case3a, 100, stacked=True)
    plt.savefig(pre+'.png')
    plt.clf()
    htmlmain.write('<tr><th>%s</th><td>%s</td><td>%s (%.f%%)</td><td>%s</td><td>%s</td><td><a href="%s.png">Open</a></td></tr>'%(d, cc_total, cc_less400, 100*cc_less400/float(cc_total), cc_more400, cc_diffchr, pre[len(workpath):]))
htmlmain.write('''
</table>

<h2>Funsion reads analysis</h2>
''')
for d, f1, f2 in dataset:
    print '> Fusion', d
    bam = workpath+'%s.%s.bowtie.fix.bam'%(d, genome)
    if not os.path.exists(bam+'.realign'):
        os.system('python realign_reads_par.py '+bam)
    if not os.path.exists(bam+'.realign.addition.png'):
        os.system('python detect_fusion_align.py '+bam+'.realign')
    htmlmain.write('<h3>%s</h3><img width=500bp src=%s.realign.addition.png></img><img width=500bp src=%s.realign.findC.png></img>'%(d, bam[len(workpath):], bam[len(workpath):]))
htmlmain.write('''

<h2> Map by Tophat fusion</h2>
<table>
<tr><th>Name</th><th>Command</th><th>Log</th></tr>
''')
for d, f1, f2 in dataset:
    print '> Tophat fusion', d, 'to', genome
    out = 'tophat_%s'%d
    cmd = 'tophat2 -o %s%s -p %s --fusion-search --keep-fasta-order --bowtie1 --no-coverage-search -r 0 --mate-std-dev 80 --max-intron-length 100000 --fusion-min-dist 100000 --fusion-anchor-length 13 /data/BOWTIE_INDEXES/%s '%(workpath, out, thread, genome) + '%s %s'%(f1, f2)
    if not os.path.exists(workpath+out+'/accepted_hits.bam'):
        os.system(cmd)
    htmlmain.write('<tr><th>%s</th><td>%s</td><td><a href="%s/align_summary.txt">Open</a></td></tr>'%(d, cmd, out))

#os.chdir(workpath)
#cmd2 = 'tophat-fusion-post -o tophatfusion_out -p %s --num-fusion-reads 1 --num-fusion-pairs 1 --num-fusion-both 1 /data/BOWTIE_INDEXES/%s '%(thread, genome)
#os.system(cmd2)
#os.chdir('..')

htmlmain.write('''
</table>

<h1>End of Report</h1>
''')

htmlmain.close()
os._exit(0)

htmlmain.write('''
<h2>Read counts by transcripts<h2>
'''%(genome, server['address'], genome))

refseq2gene = {}
if os.path.exists('refseq2gene_mm9.txt'):
    infile = open('refseq2gene_mm9.txt', 'r')
    for line in infile:
       a,b = line.split('\t') 
       refseq2gene[a] = b.strip()
    infile.close()
    print 'Get', len(refseq2gene), 'refseq ids mapped to gene names'

vals = {}
for i in xrange(len(dataset)):
    d, f1, f2 = dataset[i]
    val = {}
    temp = open(workpath+'feature_counts_brief.txt', 'r')
    temp.readline() ## skip header
    for line in temp:
        ele = line.split('\t')
        v = int(ele[i+1])
        if v >0:
            val[refseq2gene.get(ele[0], 'Unknown')+'\t'+ele[0]] = v
    temp.close()
    vals[d] = val

output = open(workpath+'feature_counts_gene.txt', 'w')
genes = val.keys()
genes.sort()
output.write('Gene\tRefseq\t%s\n'%'\t'.join([d for d,f1,f2 in dataset]))
for gene in genes:
    output.write(gene)
    for d,f1,f2 in dataset:
        val = vals[d]
        if gene in val:
            output.write('\t%s'%val[gene])
        else:
            output.write('\t0')
    output.write('\n')
output.close()

htmlmain.write('''
<table>
<tr><td>%s</td></tr>
<tr><td>%s</td></tr>
<tr><td>Output file is <a href="feature_counts_gene.txt">feature_counts_gene.txt</a></td></tr>
</table>
'''%(cmd1, cmd2))

htmlmain.write('''
<h2>Overlap of expressed genes</h2>
<table>
''')
for d1, f1, f2 in dataset:
    htmlmain.write('<tr><th>'+d1+'</th>')
    v1 = vals[d1]
    htmlmain.write('<td>%d expressed</td>'%len(v1))
    for d2, f1, f2 in dataset:
        v2 = vals[d2]
        g1 = set([g for g in v1 if v1[g]>0])
        g2 = set([g for g in v2 if v2[g]>0])
        htmlmain.write('<td>%.f%% in %s</td>'%(100*len(g1&g2)/float(len(g1)), d2))
    htmlmain.write('</tr>')

htmlmain.write('''
</table>

<h1>End of Report</h1>
''')
htmlmain.close()
