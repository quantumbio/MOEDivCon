#!/usr/bin/env python
#
# This script performs two major functions:
#   1) Submit all of the PBS jobs in newdir and wait for them to complete.
#   2) Compare the output for correctness.
#
# Created on 2010-12-07 by Scott Woods <scott@westarete.com> for QuantumBio Inc.

import os.path, re, subprocess, shlex, string, sys, time
from glob import glob
from pprint import pprint
from sets import Set

if len(sys.argv) != 3:
	print "Usage: ./dbdiff.py basedir newdir"
	sys.exit(1)
basedir = sys.argv[1]
newdir  = sys.argv[2]

# Submit each job.
submitted_job_ids = set([])
for scriptfile in glob(newdir + '/*.pbs'):
	job_id = string.strip(subprocess.Popen(['qsub', scriptfile], stdout=subprocess.PIPE).communicate()[0])
	print 'Submitted job ' + job_id
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
print 'All jobs are complete'

# A function to run a command and check the return value.
def run(cmd):
	print "Running command: " + cmd
	retval = subprocess.call(shlex.split(cmd))
	if retval != 0:
		raise Exception("command returned " + retval)

# Check the results for correctness.
for basefile in glob(basedir + '/*.mdb'):
        filename = os.path.basename(basefile)
        if re.search("pwd", filename):
		testpwd = 1
	else: 
		testpwd = 0
	newfile = newdir + '/' + filename
	if not os.path.exists(newfile):
		raise Exception('Fatal: "' + newfile + '" does not exist')
	if testpwd:
		run('moebatch -exec "run [\'pwdbatch.svl\', [\'' + newfile + '\', \'QMScore\', 15]]" -exit')
	run('moebatch -exec "run [\'qbmdbdiff.svl\', [\'qbmdbdiff-error.log\', \'' + basefile + '\', \'' + newfile + '\', \'all\', 0.5, ' + testpwd + ']]" -exit')

