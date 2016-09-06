import os
import subprocess

path = 'data/'
workpath = 'work/'
genome = 'rn5'
threads = 10

dataset = []
counts = {}
for f in os.listdir(path):
    for fq in os.listdir(path+f):
        if fq.endswith('fastq.gz'):
            dataset.append([f, path+f+'/'+fq])
dataset.sort()

for d in dataset:
    print '> Fastqc', d[0]
    out = workpath+d[-1].split('/')[-1].replace('.fastq.gz','')+'_fastqc'
    if not os.path.exists(out):
        os.system('fastqc '+d[-1]+' -t '+str(threads)+' -o '+workpath)
    output = subprocess.Popen('cat '+out+"/fastqc_data.txt | grep ^Total | awk '{print $3}'", shell=True, stdout=subprocess.PIPE).communicate()[0]
    counts[d[0]] = {'Raw reads':int(output)}

for d in dataset:
    print '> Clean', d[0]
    out = workpath+d[0]+'.trim.fastq'
    if not os.path.exists(out):
        os.system('cutadapt -a file:TruSeq2-SE.fa -a "A{100}" -m 25 '+d[-1]+' -o '+out)
    d.append(out)
    out = workpath+d[0]+'.trim.filt.fasta'
    if not os.path.exists(out):
        os.system('./CIMS/fastq_filter.pl -v -f mean:0-24:20 -of fasta '+d[-1]+' '+out)
    d.append(out)
    out = workpath+d[0]+'.trim.filt.rmdup.fasta'
    if not os.path.exists(out):
        os.system('./CIMS/fasta2collapse.pl -v '+d[-1]+' '+out)
    d.append(out)
    out = workpath+d[0]+'.trim.filt.rmdup.rmbar.fasta'
    if not os.path.exists(out):
        os.system('./CIMS/stripBarcode.pl -len 4 -format fasta '+d[-1]+' '+out) ## default -len 5
    d.append(out)

for d in dataset:
    print '> Map', d[0], 'to', genome
    out = workpath+'%s.%s.novoalign'%(d[0], genome)
    if not os.path.exists(out):
        os.system('./novoalign_pool.py -t 85 -d ref/'+genome+'.clean.idx -f '+d[-1]+' -F FA -l 25 -s 1 -r None > '+out)
    d.append(out)

for d in dataset:
    print '> Parse', d[0], 'to', genome
    out = d[-1]+'.mutation.bed'
    out2 = d[-1]+'.tag.bed'
    if not os.path.exists(out) or not os.path.exists(out2):
        os.system('./CIMS/novoalign2bed.pl -v --mismatch-file %s %s %s'%(out, d[-1], out2))
    out3 = d[-1]+'.tag.uniq.bed'
    if not os.path.exists(out3):
        os.system('./CIMS/tag2collapse.pl -v --random-barcode -EM 30 --seq-error-model alignment --weight-in-name --keep-max-score --keep-tag-name '+out2+' '+out3)
    out4 = d[-1]+'.tag.uniq.bedgraph'
    if not os.path.exists(out4):
        os.system('./CIMS/tag2profile.pl -v -ss -exact -of bedgraph -n "'+d[0]+' Uniq Tags" '+out3+' '+out4)
    out4_bw = d[-1]+'.tag.uniq.bw'
    if not os.path.exists(out4_bw):
        os.system('cat '+out3+' | sort -k1,1 -k2,2n | bedtools genomecov -bg -split -g ref/'+genome+'.chrom.sizes -i - > '+out3+'.bdg')
        os.system('bedGraphToBigWig '+out3+'.bdg ref/'+genome+'.chrom.sizes '+out4_bw)
    out5 = d[-1]+'.tag.uniq.cluster.bed'
    if not os.path.exists(out5):
        os.system('./CIMS/tag2cluster.pl -v -s -maxgap "-1" '+out3+' '+out5)
    out6 = d[-1]+'.tag.uniq.cluster.PH.bed'
    if not os.path.exists(out6):
        os.system('./CIMS/extractPeak.pl -s -v %s %s %s'%(out5, out4, out6))

for d in dataset:
    print '> CIMS', d[0], 'to', genome
    mut = d[-1]+'.mutation.bed'
    tag = d[-1]+'.tag.uniq.bed'
    out = d[-1]+'.tag.uniq.mutation.bed'
    if not os.path.exists(out):
        os.system('./CIMS/joinWrapper.py %s %s 4 4 N %s'%(mut, tag, out))
    out_del = d[-1]+'.tag.uniq.mutation_del.bed'
    if not os.path.exists(out_del):
        os.system("cat "+out+" | awk '{if($9=="+'"-"'+") {print $0}}' | cut -f 1-6 > "+out_del)
    out_del_bdg = d[-1]+'.tag.uniq.mutation_del.bedgraph'
    if not os.path.exists(out_del_bdg):
        os.system('./CIMS/tag2profile.pl -v -ss -exact -of bedgraph -n "'+d[0]+' Deletions" '+out_del+' '+out_del_bdg)
    out_del_bw = d[-1]+'.tag.uniq.mutation_del.bw'
    if not os.path.exists(out_del_bw):
        os.system('cat '+out_del+' | sort -k1,1 -k2,2n | bedtools genomecov -bg -split -g ref/'+genome+'.chrom.sizes -i - > '+out_del+'.bdg')
        os.system('bedGraphToBigWig '+out_del+'.bdg ref/'+genome+'.chrom.sizes '+out_del_bw)
    CIMS = d[-1]+'tag.uniq.mutation_del.CIMS.txt'
    cache = d[-1]+'tag.uniq.mutation_del.cache'
    if not os.path.exists(CIMS):
        os.system('./CIMS/CIMS.pl -v -n 5 -p -c %s --keep-cache %s %s %s'%(cache, tag, out_del, CIMS))
    stat = d[-1]+'tag.uniq.mutation_del.pos.dist'
    if not os.path.exists(stat):
        os.system('sort -n '+cache+'/mutation.pos.txt | uniq -c > '+stat)
    fdr = d[-1]+'tag.uniq.mutation_del.CIMS.s30.txt'
    if not os.path.exists(fdr):
        os.system("awk '{if($9<=0.001) {print $0}}' "+CIMS+" | sort -k 9,9n -k 8,8nr -k 7,7n > "+fdr)
    bed = d[-1]+'tag.uniq.mutation_del.CIMS.s30.41nt.bed'
    if not os.path.exists(bed):
        os.system("awk '"+'{print $1"\t"$2-20"\t"$3+20"\t"$4"\t"$5"\t"$6}'+"' "+fdr+" > "+bed)
    d.append(bed)
os.system('rm -rf tag2*')

for d in dataset:
    print '> MEME-CHIP', d[0]
    bed = d[-1]
    out = workpath+d[0]+'_meme'
    if not os.path.exists(out):
        os.system('bedtools getfasta -fi ref/'+genome+'.fa -s -bed '+bed+' -s -fo '+bed+'.fa')
        os.system('meme-chip '+bed+'.fa -meme-p '+str(threads)+' -ccut 100 -norand -order 1 -meme-mod oops meme-minw 4 -meme-maxw 10 -norc -old-clustering -o '+out)
        
