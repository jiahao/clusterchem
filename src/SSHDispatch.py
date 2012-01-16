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
    from SGE import SGE
    import os, subprocess

    slots = []
    for node in SGE().get_queue_instance_status():
        num_slots = max(0, node['maxslots']-int(node['load']))
        for _ in range(num_slots):
            slots.append(node['name'])

    #TODO HANDLE THE CASE WHERE WE ALREADY RAN THIS SCRIPT AND THERE MIGHT BE
    #OFF-QUEUE JOBS ALREADY RUNNING ON NODES

    num_slots = len(slots)
    jobsleft = [job for job in SGE().getuserjobs() if job.state == 'qw']
    for idx, job in enumerate(jobsleft[:num_slots]):
        print slots[idx], job.id

        procid = str(job.id)
        scriptfilename = 'exec-'+procid+'.sh'
        f = open(scriptfilename, 'w')
        #Serialize all env variables
        for var, val in os.environ.items():
            if 'QC' in var or 'PATH' in var:
                f.write("export %s=%s\n" % (var, val))

        jobdata = SGE().get_job_data(job.id)

        f.write("cd "+jobdata["sge_o_workdir"]+'\n')
        #TODO HANDLE ARRAY JOB
        #export SGE_TASK_ID=%s
        f.write("""
        set +e
/usr/bin/lamboot
        """)
        f.write(jobdata["script_file"]+' '+jobdata["job_args"].replace(","," ")+'\n')
        f.close()
        os.chmod(scriptfilename, stat.S_IEXEC | stat.S_IWRITE | stat.S_IREAD)
        
        pid = subprocess.Popen(['ssh', slots[idx], '/usr/bin/nohup',
            os.path.abspath(scriptfilename)],
            stdout=subprocess.PIPE, stderr=subprocess.PIPE,
            stdin=subprocess.PIPE)

        print ('DONE: ssh '+slots[idx]+' '+os.path.abspath(scriptfilename))
        
    exit()

    #HERE IS STUFF FOR ARRAY JOBS
    #TODO INTEGRATE
    #Scan through scripts already generated
    scriptfilenameprefix = ''
    taskinfofile = 'missing.txt'
    
    batch = []
    for line in open(taskinfofile):
        procid = line.split()[0]
        if procid in subjobs:
            scriptfilename = scriptfilenameprefix+'exec-'+procid+'.sh'
            Handle(procid)
            
