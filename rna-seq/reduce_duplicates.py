'''
Remove duplicated reads by checking the starting positions
Need to have a sorted SAM file
Example:

samtools view -h in.bam | python reduce_duplicates.py 4 | samtools -bS - -o out.bam
'''
import sys

max_count = int(sys.argv[1])

last = ''
count = 0
for line in sys.stdin:
    if line.startswith('@'): ## header
        sys.stdout.write(line)
        continue
    ele = line.split()
    now = ele[2]+':'+ele[3]
    flag = bin(int(ele[1]))[3:]
    if len(flag) >= 3 and flag[-3] == '1': ## unmapped
        continue
    if now == last:
        count += 1
        if count <= max_count:
            sys.stdout.write(line)
    else:
        last = now
        count = 1
        sys.stdout.write(line)
