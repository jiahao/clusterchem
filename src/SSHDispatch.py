#!/usr/bin/env python

"""
Manual jobs dispatcher via SSH - PURE EVIL.

This script handles the unfortunate case where the hard coded queue quota
results in lots of idle machines when no one else is running jobs.

Tell this script what machine and which array task IDs to initiate. 
"""

import os
import stat
from multiprocessing import Pool

def Handle(procid):
    print 'Running', procid
    scriptfilename = scriptfilenameprefix+'exec-'+procid+'.sh'
    f = open(scriptfilename, 'w')
    #Serialize all env variables
    for var, val in os.environ.items():
        if 'QC' in var or 'PATH' in var:
           f.write("export %s=%s\n" % (var, val))
    f.write("""
cd %s
export SGE_TASK_ID=%s
set +e
/usr/bin/lamboot
python ~/rmt/bin/JobHandler.py --taskfile missing.txt
         """ % (os.getcwd(), str(procid)))
    f.close()
    os.chmod(scriptfilename, stat.S_IEXEC | stat.S_IWRITE | stat.S_IREAD)
    print ('ssh '+machine+' '+os.path.abspath(scriptfilename))
    import subprocess
    pid = subprocess.Popen(['ssh', machine, '/usr/bin/nohup', os.path.abspath(scriptfilename)], stdout=subprocess.PIPE, stderr=subprocess.PIPE, stdin=subprocess.PIPE)
    print ('DONE: ssh '+machine+' '+os.path.abspath(scriptfilename))

if __name__ == '__main__':
    import sys
    machine = sys.argv[1]
    
    subjobs = sys.argv[2:]
    
    #Scan through scripts already generated
    scriptfilenameprefix = ''
    taskinfofile = 'missing.txt'
    
    batch = []
    for line in open(taskinfofile):
        procid = line.split()[0]
        if procid in subjobs:
            scriptfilename = scriptfilenameprefix+'exec-'+procid+'.sh'
            Handle(procid)
            
