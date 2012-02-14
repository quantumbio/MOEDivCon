#!/usr/bin/env python
#
#   Compare the output for correctness.
#
# Created on 2010-12-07 by Scott Woods <scott@westarete.com> for QuantumBio Inc.

import os.path, re, subprocess, shlex, sys
from glob import glob

if len(sys.argv) != 3:
    print "Usage: ./dbdiff.py basedir newdir"
    sys.exit(1)
basedir = sys.argv[1]
newdir  = sys.argv[2]

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
        testpwd = '1'
    else: 
        testpwd = '0'
    newfile = newdir + '/' + filename
    if not os.path.exists(newfile):
        raise Exception('Fatal: "' + newfile + '" does not exist')
    if testpwd:
        run('qbmoebatch -exec "run [\'pwdbatch.svl\', [\'' + newfile + '\', \'QMScore\', 15]]" -exit')
    run('qbmoebatch -exec "run [\'qbmdbdiff.svl\', [\'qbmdbdiff-error.log\', \'' + basefile + '\', \'' + newfile + '\', \'all\', 0.5, ' + testpwd + ']]" -exit')
output_file = "qbmdbdiff-error.log"
if not os.path.exists(output_file):
    print "Output file does not exist. Something went wrong with the comparison."
if os.path.getsize(output_file) == 0:
    print "Output file is empty. Databases are the same."
else:
    print "Databases differ. See https://ci.quantumbioinc.com/job/MOEDivconIntegrationTests/ws/MOEDivcon/functests/qbmdbdiff-error.log for details."
    print "It's never there.  Instead it is at /home/hudson/shared_workspace/$JOB_NAME/MOEDivcon/functests.  Studying how this after failure can be put back into workspace. Starts in svl"
    sys.exit(1)
