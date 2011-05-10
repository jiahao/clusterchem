#!/usr/bin/env python

"""
Jiahao Chen <jiahao@mit.edu> 2011-02-18

This script builds heavily upon the QM-MM wrapper scripts of Tim Kowalczyk, Lee-Ping Wang and Shane Yost.

.. TODO:: UNTESTED.
"""

import glob
import os
import sys
import signal


def PrintHelp():
    print "usage: "+sys.argv[0]+" [CHARMM input file name] [additional CHARMM input keywords]"



def copywild(src, dest): #NOT recursive
    import shutil
    for filename in glob.glob(src):
        if os.path.isfile(filename):
            try:
                shutil.copy(filename, os.path.join(src, dest))
            except IOError, e:
                print 'Warning, skipped',filename
                print 'Error was:',e




def UnboundCHARMMVariables(InputFile, OtherInputs = ''):
    """Checks a CHARMM input file for unbound variables (i.e. variables
that are used but not defined.

OtherInputs is an optional input that is additional command line input.

"""

    namespace = [x.split('=')[0] for x in OtherInputs.upper().split()]
    for l in open(InputFile):
        t = l.upper().split()
        if len(t) > 0 and t[0] == 'SET':
            var = t[1].split('=')[0]
            if var not in namespace:
                namespace.append(var)

    unbound = []
    for l in open(InputFile):
        for t in l.upper().split():
            if '@' in t:
                var = t[t.index('@')+1:].replace(')','').replace('"','')
                if var not in namespace and var not in unbound:
                    unbound.append(var)

    return unbound



def RunCHARMM(InputFile, OutputFile, OtherInputs = '', CHARMMBin = 'charmm', Overwrite = True, PollInt = 60):


    if not Overwrite and os.path.exists(OutputFile):
        print 'Output file', OutputFile, 'already exists. Abort.'
        return

    
    UnboundVars = UnboundCHARMMVariables(InputFile, OtherInputs)
    if len(UnboundVars) > 0:
        assert False, 'ERROR: Unbound variables detected: '+' '.join(UnboundVars)
    
    from subprocess import Popen, PIPE, STDOUT

    os.setpgrp() #create new process group, become its leader

    try:
        p = Popen(CHARMMBin+' '+OtherInputs, stdin=PIPE, stdout=PIPE, stderr=STDOUT, bufsize=1, shell=True, close_fds=True)
        #p.stdin.write('\n'.join([open(InputFile).read(), OtherInputs]))
        p.stdin.write(open(InputFile).read())
        p.stdin.close()
    
        o = open(OutputFile, 'w', 1)
        o.write(p.stdout.read())
        p.stdout.close()
        o.close()
    except KeyboardInterrupt:
        #Kill ALL children
        os.killpg(0, signal.SIGKILL)
        exit()
    



if __name__ == '__main__':
    if len(sys.argv) == 1:
        PrintHelp()
        exit()

    #Make a scratch directory in user scratch space
    ScratchDir = os.path.join('/scratch', os.environ['USER'], str(os.getpid())+'.qmmmpol')
    if not os.path.exists(ScratchDir):
        os.makedirs(ScratchDir)
    #Stage scratch directory
    copywild('*', ScratchDir)
    #Save scratch dir in environment
    os.environ['QMMMPOL_SCR'] = ScratchDir

    CHARMMInput = sys.argv[1]
    CHARMMFileRoot = '.'.join(CHARMMInput.split('.')[:-1])
    CHARMMOutput= CHARMMFileRoot + '.out'

    #Save current PID
    os.environ['QMMMPOL_PID'] = str(os.getpid())

    #Save working directory
    WorkingDir = os.getcwd()
    os.environ['QMMMPOL_OWD'] = WorkingDir


    OtherInputs = ' '.join(sys.argv[2:])

    #Run CHARMM
    RunCHARMM(CHARMMInput, os.path.join(WorkingDir, CHARMMOutput), OtherInputs)
