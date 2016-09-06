import os

path = 'data/'
workpath = 'work/'
genome = 'hg19'
threads = 25
max_trim = 6
RBP = 'nPTB'

dataset = []
for f in os.listdir(path):
    if path.startswith('test'):
        name = f.replace('.fastq.gz','')
        dataset.append((name, path+f))
        continue
    if f.find(RBP)==0 and f.find('1')>0 and f.find('ko')<0 and f.endswith('fastq.gz'):
        name = f.replace('.fastq.gz','')
        dataset.append((name, path+f))
dataset.sort()

for d,f in dataset:
    print '> Check', d, 'to', genome
    out = workpath+'%s.fastx.log'%d
    if not os.path.exists(out):
        os.system('zcat ' + f + ' | fastx_quality_stats >' + out)

for d,f in dataset:
    print '> Map', d, 'to', genome
    out = workpath+'%s.%s.bowtie.bam'%(d, genome)
    if not os.path.exists(out):
        for i in xrange(0, max_trim+1):
            if i == 0:
                cmd = 'zcat '+f+'|'
            else:
                cmd = 'cat '+workpath+'temp_bowtie_output%s.fq'%(i-1)+'|'
            filt = ' -F 4 ' ## only for mapped read pairs
            os.system(cmd+'bowtie /data/BOWTIE_INDEXES/'+genome+' -p '+str(threads)+' --trim5 '+str(i)+' -l25 -n2 -k11 -m10 -e200 --best --strata --sam --phred33-quals - --un '+workpath+'temp_bowtie_output'+str(i)+'.fq | samtools view -ubS'+filt+'- -o '+workpath+'temp_bowtie_%s.bam'%i)
        os.system('samtools merge -c '+out+' '+workpath+'temp_bowtie_*.bam')
        os.system('samtools flagstat '+out+' >'+out+'.log')
        os.system('rm '+workpath+'temp_bowtie*')

for d,f in dataset:
    print '> Sort and rmdup', d, 'to', genome
    pre = workpath+'%s.%s.bowtie'%(d, genome)
    if not os.path.exists('%s.rmdup.bam.bai'%pre):
        os.system('samtools sort -@ %s %s.bam -o %s.srt.bam'%(threads, pre, pre)) ## sort bam
        os.system('samtools view -h %s.srt.bam | python reduce_duplicates.py | samtools view -bS -o %s.rmdup.bam'%(pre, pre)) ## remove dup
        os.system('samtools index %s.rmdup.bam'%pre)

sites = []
for d,f in dataset:
    print '> MiClip', d, 'to', genome
    inp = workpath+'%s.%s.bowtie.rmdup.bam'%(d, genome)
#    out = inp+'.enrich'
    out = inp+'.miclip'
    if not os.path.exists(out):
        os.system('samtools view %s -o %s.sam'%(inp, inp))
        sfile = open(out+'.r', 'w')
        sfile.write('''
library("MiClip")
mc = MiClip("%s.sam")
rd = MiClip.read(mc)
en = MiClip.enriched(rd)
bd = MiClip.binding(en)
#tag = aggregate(tag ~ region_id, data=bd$sites, max)
#mut = aggregate(mutant ~ region_id, data=bd$sites, max)
#out = merge(bd$cluster, tag, by="region_id")
#out = merge(out, mut, by="region_id")
out = bd$sites
write.table(out, "%s", sep="\\t")
        '''%(inp, out))
        sfile.close()
        os.system('/usr/bin/Rscript %s.r'%out)
    infile = open(out, 'r')
    infile.readline()
    for line in infile:
        ele = line.split('\t')
#        if ele[8] == 'TRUE':
#            sites.append((ele[4].replace('"',''), int(float(ele[5]))-50, int(float(ele[5]))+50, ele[3].replace('"',''), '%s|%s'%(d, ele[0].replace('"','')), int(ele[7])))
        if ele[6] == 'TRUE' and ele[7] == 'TRUE' and int(ele[8]) >= 10:
            sites.append((ele[2].replace('"',''), int(ele[4]), int(ele[5]), ele[3].replace('"',''),'%s|%s'%(d,ele[1]), int(ele[8])))
    infile.close()
sites.sort()
print 'We have', len(sites), 'binding sites.'
output = open(workpath+'miclip_sites_'+RBP+'.bed', 'w')
for ch,p1,p2,st,nm,sc in sites:
    output.write('\t'.join([ch,str(p1),str(p2),nm,str(sc),st])+'\n')
output.close()

top_sites = sorted(sites, key=lambda x: x[5], reverse=True)
output = open(workpath+'miclip_top_sites_'+RBP+'.bed', 'w')
for ch,p1,p2,st,nm,sc in top_sites[:1000]:
    output.write('\t'.join([ch,str(p1),str(p2),nm,str(sc),st])+'\n')
output.close()

#os.system('bedtools getfasta -fi /data/BOWTIE_INDEXES/'+genome+'.fa -s -bed '+workpath+'miclip_top_sites_'+RBP+'.bed -s -name -fo '+workpath+'miclip_top_sites_'+RBP+'.fa')
#os.system('rm -rf '+workpath+'motif_'+RBP)
#os.system('meme-chip '+workpath+'miclip_top_sites_'+RBP+'.fa -meme-p '+str(threads)+' -nmeme 1000 -ccut 100 -norand -order 1 -meme-mod oops -meme-nmotifs 5 -meme-minw 6 -meme-maxw 8 -norc -old-clustering -o '+workpath+'motif_'+RBP)
#os.system('findMotifsGenome.pl '+workpath+'miclip_top_sites_'+RBP+'.bed '+genome+' '+workpath+'motif_'+RBP+' -bits -redundant 0.9 -rna -size 20 -len 6')

os.system('bedtools intersect -s -u -a refseq_'+genome+'.bed -b '+workpath+'miclip_sites_'+RBP+'.bed > '+workpath+'miclip_igene.bed')
os.system('bedtools intersect -s -v -a refseq_'+genome+'.bed -b '+workpath+'miclip_sites_'+RBP+'.bed > '+workpath+'miclip_ogene.bed')

def filter_bed(inbed, outbed, min_len=2000):
    os.system('bedtools intersect -wa -wb -s -a '+inbed+' -b '+inbed+' > '+inbed+'tools.bed')
    inf = open(inbed+'tools.bed', 'r')
    out = open(outbed, 'w')
    gene = {}
    link = {}
    for line in inf:
        ele = line.split()
        if ele[3] == ele[9]:
            if ele[3] in gene:
                e1 = gene[ele[3]]
                l1 = abs(int(e1[1])-int(e1[2]))
                length = abs(int(ele[1])-int(ele[2]))
                if l1 > length:
                    continue
            gene[ele[3]] = tuple(ele[:6])
        else:
            g1 = ele[3]
            g2 = ele[9]
            gp = link.get(g1, [])
            gp.append(g2)
            link[g1] = gp
    mgene = []
    for g1 in gene:
        if g1 not in link or len(link[g1]) == 1:
            mgene.append(g1)
            continue
        e1 = gene[g1]
        l1 = abs(int(e1[1])-int(e1[2]))
        for g2 in link[g1]:
            e2 = gene[g2]
            l2 = abs(int(e2[1])-int(e2[2]))
            if l2 > l1:
                break
        if l1 >= l2 and g1 < g2:
            mgene.append(g1)
    mgene.sort()
    cc = 0
    for g in mgene:
        ele = gene[g]
        length = abs(int(ele[1])-int(ele[2]))
        if length < min_len:
            continue
        out.write('\t'.join(ele)+'\n')
        cc += 1
    print len(gene), len(mgene), cc
    inf.close()
    out.close()

filter_bed(workpath+'miclip_igene.bed', workpath+'miclip_igene_filt.bed')
filter_bed(workpath+'miclip_ogene.bed', workpath+'miclip_ogene_filt.bed')

#os.system('python bed2uniq.py miclip_igene.bed 1 miclip_igene_uniq.bed')
#os.system('python bed2uniq.py miclip_ogene.bed 1 miclip_ogene_uniq.bed')

#os.system('python length_filter.py miclip_igene_uniq.bed -f bed -t 2000 --outfile miclip_igene_uniq_gt2k.bed')
#os.system('python length_filter.py miclip_ogene_uniq.bed -f bed -t 2000 --outfile miclip_ogene_uniq_gt2k.bed')

#os.system('cp -f miclip_igene_uniq_gt2k.bed '+workpath+'miclip_igene.bed')
#os.system('cp -f miclip_ogene_uniq_gt2k.bed '+workpath+'miclip_ogene.bed')

os.chdir(workpath)

mlist = []
for d,f in dataset:
    mlist.append('%s.%s.bowtie.rmdup.bam'%(d, genome))
os.system('samtools merge -f -c %s.bam %s'%(RBP, ' '.join(mlist)))
os.system('samtools index %s.bam'%RBP)

with open('ngs.plot.config.txt', 'w') as conf:
    conf.write(RBP+'.bam\t"miclip_igene_filt.bed"\t"Targets"\n')
    conf.write(RBP+'.bam\t"miclip_ogene_filt.bed"\t"Non-targets"\n')

os.system("/usr/bin/Rscript /data/pipelines/softwares/bin/ngs.plot.r -AL spline -G hg19 -R bed -SS same -C ngs.plot.config.txt -O plot_%s_same -T %s_same -L 2000 -FL 40 -P 0 -MQ 10 -RB 0.01"%(RBP, RBP))
os.system("/usr/bin/Rscript /data/pipelines/softwares/bin/ngs.plot.r -AL spline -G hg19 -R bed -SS opposite -C ngs.plot.config.txt -O plot_%s_oppo -T %s_oppo -L 2000 -FL 40 -P 0 -MQ 10 -RB 0.01"%(RBP, RBP))
os.system("../replot.r prof -I plot_%s_oppo.zip -O %s_oppo.prof"%(RBP, RBP))

