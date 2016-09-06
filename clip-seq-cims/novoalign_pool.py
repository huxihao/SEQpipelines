#!/usr/bin/env python

import os, sys, shutil
from multiprocessing import Pool

def run(cmd):
    os.system(cmd+'.tmp')
    F = cmd.split()[-1]
    os.rename(F+'.tmp', F)

def start(dfile, path, groups):
    cc = 0
    inp = open(dfile, 'r')
    for line in inp:
        cc += 1
    inp.close()

    step = max(10000, (cc/groups/4)*4)
    outs = []
    inp = open(dfile, 'r')
    cc = 0
    gg = 1
    ofile = path+str(gg)+'.fa'
    outf = open(ofile, 'w')
    for line in inp:
        outf.write(line)
        cc += 1 
        if cc == step:
            outf.close()
            outs.append(ofile)
            gg += 1
            ofile = path+str(gg)+'.fa'
            outf = open(ofile, 'w')
            cc = 0
    outf.close()
    if cc > 0:
        outs.append(ofile)
    inp.close()
    return outs

if __name__ == '__main__':
    arg = sys.argv
    threads = 25
    groups = threads*20
    indx = 0
    for i in xrange(len(arg)):
        if arg[i] == '-f':
            indx = i+1

    tmp = arg[indx]+'.tmp/'
    if not os.path.exists(tmp):
        os.mkdir(tmp)

    inputs = start(arg[indx], tmp, groups)

    cmds = []
    for F in inputs:
        if not os.path.exists(F+'.out'):
            cmd = arg
            cmd[indx] = F
            cmds.append(' '.join(['novoalign']+cmd[1:]+['>', F+'.out']))

    p = Pool(threads)
    p.map(run, cmds, chunksize=1)

    os.system('cat '+' '.join([F+'.out' for F in inputs]))
#    shutil.rmtree(tmp)

