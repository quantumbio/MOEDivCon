#!/usr/bin/env python
#
#   Submit all PBS jobs in current directory and wait for them to complete.
#
# Created on 2010-12-07 by Scott Woods <scott@westarete.com> for QuantumBio Inc.

import subprocess, string, sys, time
from glob import glob

# Submit each job.
submitted_job_ids = set([])
for scriptfile in glob('*.pbs'):
    job_id = string.strip(subprocess.Popen(['qsub', scriptfile], stdout=subprocess.PIPE).communicate()[0])
    print 'Submitted job ' + scriptfile + ' with PBS id ' + job_id
    submitted_job_ids.add(job_id)

# Wait until each job completes.
started_at = time.clock()
completed_job_ids = set([])
while submitted_job_ids != completed_job_ids:
    time.sleep(5)
    for job_id in submitted_job_ids - completed_job_ids:
        retval = subprocess.call(['qstat', job_id], stdout=open('/dev/null', 'w'))
        if retval > 0:
            print 'Job ' + job_id + ' is no longer running'
            completed_job_ids.add(job_id)
    hour_limit = 12
    if submitted_job_ids != completed_job_ids and time.clock() - started_at > 60*60*hour_limit:
        print "Fatal: Waited " + hour_limit + " hours for all jobs to complete."
        print "Here are the ones that are still running:"
        for job_id in submitted_job_ids - completed_job_ids:
            print "    " + job_id
        sys.exit(1)
print 'All jobs are complete'
