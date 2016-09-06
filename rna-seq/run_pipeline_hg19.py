import os
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--datapath', type=str, default='YU20160627/', help='Data path')
parser.add_argument('--workpath', type=str, default='work_strand/', help='Work path')
parser.add_argument('--genome', type=str, default='hg19', help='Reference genome')
parser.add_argument('--threads', type=int, default=15, help='Number of threads')
parser.add_argument('--server', type=str, default='../server.txt', help='File for saving the server information')

args = parser.parse_args()
print args

genome = args.genome
threads = args.threads
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
    if len(f) == 2:
        dataset.append((p, datapath+p+'/'+f[0], datapath+p+'/'+f[1]))
dataset.sort()

htmlmain.write('''
<h1>Report for RNA-seq</h1>
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
        os.system('fastqc -t '+str(threads)+' '+f1+' -o '+workpath)
    f2name = (f2.split('/')[-1]).split('.')[0]+'_fastqc'
    if not os.path.exists(workpath+f2name):
        os.system('fastqc -t '+str(threads)+' '+f2+' -o '+workpath)
    htmlmain.write('<tr><th>%s</th><td><a href="%s">%s</a></td><td><a href="%s/fastqc_report.html" target="_blank">Read 1</a> <a href="%s/fastqc_report.html" target="_blank">Read 2</a></td></tr>'%(d, out, out, f1name, f2name))
htmlmain.write('''
</table>        
        
<h2>Remove adapters</h2>
<table>
<tr><th>Name</th><th>Command</th><th>Fastqc output</th></tr>
 ''')
newset = []
for d, f1, f2 in dataset:
    print '> Trim', d, 'to', genome
    F1 = 'clean/'+d+'_R1.fq.gz'
    F2 = 'clean/'+d+'_R2.fq.gz'
    cmd = 'java -jar trimmomatic/trimmomatic-0.33.jar PE '+f1+' '+f2+' '+F1+' /dev/null '+F2+' /dev/null ILLUMINACLIP:trimmomatic/adapters/TruSeq3-PE-2.fa:2:30:10 LEADING:6 TRAILING:20 SLIDINGWINDOW:4:15 MINLEN:36'
    if not os.path.exists(F1) or not os.path.exists(F2):
        os.system(cmd)
    f1name = d+'_R1.fq_fastqc'
    if not os.path.exists(workpath+f1name):
        os.system('fastqc -t '+str(threads)+' '+F1+' -o '+workpath)
    f2name = d+'_R2.fq_fastqc'
    if not os.path.exists(workpath+f2name):
        os.system('fastqc -t '+str(threads)+' '+F2+' -o '+workpath)
    newset.append((d, F1, F2))
    htmlmain.write('<tr><th>%s</th><td>%s</td><td><a href="%s/fastqc_report.html" target="_blank">Read 1</a> <a href="%s/fastqc_report.html" target="_blank">Read 2</a></td></tr>'%(d, cmd, f1name, f2name))
dataset = newset

htmlmain.write('''
</table>        

<h2>Map reads by Tophat</h2>
<table>
<tr><th>Name</th><th>Command</th><th>Log</th></tr>
''')
for d, f1, f2 in dataset:
    print '> Tophat', d, 'to', genome
    out = '%s.%s.tophat.bam'%(d, genome)
    cmd = 'tophat2 -G refseq_%s.gtf --transcriptome-index=transidx --library-type fr-secondstrand -p %s %s '%(genome, threads, genome) + '%s %s'%(f1,f2)
    if not os.path.exists(workpath+out):
        os.system(cmd)
        os.system('mv tophat_out/accepted_hits.bam '+workpath+out)
        os.system('mv tophat_out/junctions.bed '+workpath+out+'.junctions.bed')
        os.system('rm -rf tophat_out')
    if not os.path.exists(workpath+out+'.log'):
        os.system('samtools flagstat '+workpath+out+'>'+workpath+out+'.log')
    htmlmain.write('<tr><th>%s</th><td>%s</td><td><a href="%s.log">Open</a></td></tr>'%(d, cmd, out))
    if not os.path.exists(workpath+out+'.bai'):
        os.system('samtools index '+workpath+out)
htmlmain.write('''
</table>

<h2>Remove duplicates<h2>
<table>
<tr><th>Name</th><th>Before</th><th>After</th></tr>
''')
#for d, f1, f2 in dataset:
#    print '> rmdup', d
#    inp = workpath+'%s.%s.tophat.bam'%(d, genome)
#    out = workpath+'%s.%s.tophat.rmdup.bam'%(d, genome)
#    if not os.path.exists(out):
#        os.system('samtools view -h %s | python reduce_duplicates.py 4 | #samtools view -bS - -o %s'%(inp, out))
#    if not os.path.exists(workpath+inp+'.log'):
#        os.system('samtools flagstat '+inp+'>'+inp+'.log')
#    with open(inp+'.log', 'r') as tmp:
#        val1 = tmp.readline().split()[0]
#    if not os.path.exists(workpath+out+'.log'):
#        os.system('samtools flagstat '+out+'>'+out+'.log')
#    with open(out+'.log', 'r') as tmp:
#        val2 = tmp.readline().split()[0]
#    htmlmain.write('<tr><th>%s</th><td>%s</td><td>%s</td></tr>'%(d, val1, val2))
htmlmain.write('''
</table>

<h2>Genome browser tracks</h2>
<p>Copy the following tracks:</p>
<table><tr><td>''')
for d, f1, f2 in dataset:
    print '> BAM to BigWig', d
    out1 = '%s.%s.tophat_pos.bw'%(d, genome)
    out2 = '%s.%s.tophat_neg.bw'%(d, genome)
    htmlmain.write('track type=bigWig db=%s name="%s Pos" visibility=2 bigDataUrl="%s/%s/rna-seq/%s"<br/>'%(genome, d, server['address'], server['home'], out1))
    htmlmain.write('track type=bigWig db=%s name="%s Neg" visibility=2 bigDataUrl="%s/%s/rna-seq/%s"<br/>'%(genome, d, server['address'], server['home'], out2))
    if not os.path.exists(workpath+out2):
        os.system("bedtools genomecov -bg -split -strand + -ibam "+workpath+"%s.%s.tophat.bam > tmp.bdg"%(d, genome))
        os.system("bedGraphToBigWig tmp.bdg %s.chrom.sizes "%genome+workpath+out1)
        os.system("bedtools genomecov -bg -split -strand - -ibam "+workpath+"%s.%s.tophat.bam > tmp.bdg"%(d, genome))
        os.system("bedGraphToBigWig tmp.bdg %s.chrom.sizes "%genome+workpath+out2)
        os.system("rm tmp.bdg")

htmlmain.write('''</td></tr></table>
<p>into <a href="http://genome.ucsc.edu/cgi-bin/hgCustom?db=%s" target="_black">the UCSC genome browser</a>, or use <a href="%s/internal/browser/cgi-bin/hgCustom?db=%s" target="_black">the local genome browser</a> (user name: xuelab; password: neibu)</p>

<h2>Read counts by transcripts<h2>
<table>
<tr><th>Name</th><th>Command</th></tr>
'''%(genome, server['address'], genome))

cmd2 = 'sed 1d %sfeature_counts | cut -f1,7- > %sfeature_counts_brief.txt'%(workpath,workpath)

refseq2gene = {}
if os.path.exists('refseq2gene_%s.txt'%genome):
    infile = open('refseq2gene_%s.txt'%genome, 'r')
    for line in infile:
       a,b = line.split('\t') 
       refseq2gene[a] = b.strip()
    infile.close()
    print 'Get', len(refseq2gene), 'refseq ids mapped to gene names'

vals = {}
for d, f1, f2 in dataset:
    bam = workpath+'%s.%s.tophat.bam'%(d, genome)
    out = workpath+d+'.counts.txt'
    cmd = 'featureCounts -t exon -g gene_id -O -M -a refseq_%s.gtf -o %s %s'%(genome, out, bam)
    htmlmain.write('<tr><th>%s</th><td>%s</td></tr>'%(d, cmd))
    if not os.path.exists(out):
        os.system(cmd)
    temp = open(out, 'r')
    temp.readline()
    print temp.readline(), ## skip header
    length = {}
    val = {}
    for line in temp:
        ele = line.split('\t')
        v = int(ele[-1])
        l = int(ele[-2])
        g = refseq2gene.get(ele[0], 'Unknown')+'\t'+ele[0]
        length[g] = l
        if v > 0:
            val[g] = v
    temp.close()
    vals[d] = val
    if 'length' not in vals:
        vals['length'] = length

output = open(workpath+'feature_counts_gene.txt', 'w')
genes = vals['length'].keys()
genes.sort()
output.write('Gene\tRefseq\tLength\t%s\n'%'\t'.join([d for d,f1,f2 in dataset]))
for gene in genes:
    output.write(gene)
    length = vals['length']
    output.write('\t%s'%length[gene])
    for d,f1,f2 in dataset:
        val = vals[d]
        if gene in val:
            output.write('\t%s'%val[gene])
        else:
            output.write('\t0')
    output.write('\n')
output.close()

htmlmain.write('''
<tr><th>Overall</th><td>Summary file is <a href="feature_counts_gene.txt">feature_counts_gene.txt</a></td></tr>
</table>
''')

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
