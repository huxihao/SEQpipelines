inp = open('mm9.fa', 'r')
out = open('mm9.clean.fa', 'w')
save = False
for line in inp:
    if line[0] == '>':
        if line.find('_') < 0 and line.find('M') < 0:
            print line,
            save = True
        else:
            save = False
    if save:
        out.write(line)
out.close()
inp.close()
