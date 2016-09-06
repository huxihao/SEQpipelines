'''
Remove duplicated reads by checking the starting positions
Need to have a sorted SAM file
Example:

samtools view -h in.bam | python reduce_duplicates.py | samtools -bS - -o out.bam

'''
import sys

max_count = 4

last = ''
head = '' ## the first 4 bases are the barcode
count = 0
for line in sys.stdin:
    if line.startswith('@'): ## header
        sys.stdout.write(line)
        continue
    ele = line.split()
    now = ele[2]+':'+ele[3]
    bar = ele[0].split('|')[1]
    assert len(bar) == 4
    ## the bar code should be NNNT
    if bar.find('N') >= 0:
        continue
    if bar[-1] != 'T':
        continue
    flag = bin(int(ele[1]))[3:]
    if len(flag) >= 3 and flag[-3] == '1': ## unmapped
        continue
    if now == last and bar == head:
        count += 1
        if count <= max_count:
            sys.stdout.write(line)
    else:
        last = now
        head = bar
        count = 1
        sys.stdout.write(line)
