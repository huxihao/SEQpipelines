import os

path = 'data/'
workpath = 'work_hg18/'
genome = 'hg18'
thread = 25
max_trim = 5
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

#new_f = []
#for d,f in dataset:
#    print '> Trim', d, 'to', genome
#    out = workpath+'%s.cutadapt.fq.gz'%d
#    if not os.path.exists(out):
#        os.system('cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -o %s %s'%(out, f))
#    new_f.append((d, out))
#dataset = new_f

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
            os.system(cmd+'bowtie /data/BOWTIE_INDEXES/'+genome+' -p '+str(thread)+' --trim5 '+str(i)+' -l25 -n2 -k11 -m10 -e200 --best --strata --sam --phred33-quals - --un '+workpath+'temp_bowtie_output'+str(i)+'.fq | samtools view -ubS'+filt+'- -o '+workpath+'temp_bowtie_%s.bam'%i)
        os.system('samtools merge -c '+out+' '+workpath+'temp_bowtie_*.bam')
        os.system('samtools flagstat '+out+' >'+out+'.log')
        os.system('rm '+workpath+'temp_bowtie*')

for d,f in dataset:
    print '> sort and rmdup', d, 'to', genome
    pre = workpath+'%s.%s.bowtie'%(d, genome)
    if not os.path.exists('%s.srt.bam'%pre):
        os.system('samtools sort -@ %s %s.bam -o %s.srt.bam'%(thread, pre, pre)) ## sort bam
    if not os.path.exists('%s.rmdup.bam.bai'%pre):
        os.system('samtools view -h %s.srt.bam | python reduce_test.py | samtools view -bS -o %s.rmdup.bam'%(pre, pre)) ## remove dup
        os.system('samtools index %s.rmdup.bam'%pre)

sites = {}
for d,f in dataset:
    print '> MiClip', d, 'to', genome,
    inp = workpath+'%s.%s.bowtie.rmdup.bam'%(d, genome)
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
write.table(bd$clusters, "%s", sep="\\t")
        '''%(inp, out))
        sfile.close()
        os.system('/usr/bin/Rscript %s.r'%out)
    infile = open(out, 'r')
    infile.readline()
    site = []
    for line in infile:
        ele = line.split('\t')
        id = ele[1]
        ch = ele[2].replace('"','')
        st = ele[3].replace('"','')
        p1 = int(ele[4])
        p2 = int(ele[5])
        enriched = (ele[6].strip().lower()=='true')
        bound = (ele[7].strip().lower()=='true')
        if enriched and bound:
            site.append((ch,p1,p2,st,'%s|%s'%(d,id)))
    site.sort()
    infile.close()
    print len(site)
    sites[d] = site     

output = open(workpath+'miclip_sites.bed', 'w')
for siite in sites:
    for ch,p1,p2,st,nm in site:
        output.write('\t'.join([ch,str(p1),str(p2),nm,'1',st])+'\n')
output.close()

os.system('bedtools intersect -s -u -a refseq_'+genome+'.bed -b '+workpath+'miclip_sites.bed > '+workpath+'miclip_igene.bed')
os.system('bedtools intersect -s -v -a refseq_'+genome+'.bed -b '+workpath+'miclip_sites.bed > '+workpath+'miclip_ogene.bed')

def filter_bed(inbed, outbed, min_len=2000):
    os.system('bedtools intersect -wa -wb -s -a '+inbed+' -b '+inbed+' > '+inbed+'tools.bed')
    inf = open(inbed+'tools.bed', 'r')
    out = open(outbed, 'w')
    known = set()
    for line in inf:
        ele = line.split()
        length = abs(int(ele[1])-int(ele[2]))
        mark = '%s:%s-%s:%s'%(ele[0], ele[1], ele[2], ele[5])
        if length < min_len or mark in known:
            continue
        known.add(mark)
#        if length < min_len or ele[3] in known or ele[9] in known:
#            known.add(ele[3]) ## self
#            known.add(ele[9]) ## overlapped
#            continue
#        known.add(ele[3]) ## self
#        known.add(ele[9]) ## overlapped
        out.write('\t'.join(ele[:6])+'\n')
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

for bam in mlist:
    os.system('bedtools genomecov -ibam '+bam+' -strand + -bg > temp.bdg')
    os.system('sort -k1,1 -k2,2n temp.bdg > temp.srt.bdg')
    os.system('bedGraphToBigWig temp.srt.bdg ../'+genome+'.chrom.sizes '+bam[:-4]+'.pos.bw')
    os.system('bedtools genomecov -ibam '+bam+' -strand - -bg > temp.bdg')
    os.system('sort -k1,1 -k2,2n temp.bdg > temp.srt.bdg')
    os.system('bedGraphToBigWig temp.srt.bdg ../'+genome+'.chrom.sizes '+bam[:-4]+'.neg.bw')
    os.system('rm temp.bdg temp.srt.bdg')

with open('ngs.plot.config.txt', 'w') as conf:
    conf.write(RBP+'.bam\t"miclip_igene_filt.bed"\t"Targets"\n')
    conf.write(RBP+'.bam\t"miclip_ogene_filt.bed"\t"Non-targets"\n')

os.system('samtools merge -f -c %s.bam %s'%(RBP, ' '.join(mlist)))
os.system('samtools index %s.bam'%RBP)

os.system("/usr/bin/Rscript /data/pipelines/softwares/bin/ngs.plot.r -AL spline -G "+genome+" -R bed -SS same -C ngs.plot.config.txt -O plot_%s_same -T %s_same -L 2000 -FL 40 -P 0 -MQ 10 -RB 0.01"%(RBP, RBP))
os.system("/usr/bin/Rscript /data/pipelines/softwares/bin/ngs.plot.r -AL spline -G "+genome+" -R bed -SS opposite -C ngs.plot.config.txt -O plot_%s_oppo -T %s_oppo -L 2000 -FL 40 -P 0 -MQ 10 -RB 0.01"%(RBP, RBP))
os.system("../replot.r prof -I plot_%s_oppo.zip -O %s_oppo.prof"%(RBP, RBP))

