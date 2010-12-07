#!/bin/bash
#

set -e  # Exit if we encounter any errors.
set -x  # Trace each command as it's executed.

if [ $# != 2 ] ; then
  echo "Usage:  dbdiff.sh basedir newdir"
  exit
fi

#  Submit new jobs to PBS queue

basefiles=`/bin/ls -1 $1/*.mdb`

for bf in $basefiles ; do
  bn=`basename $bf`
  nf="${2}/$bn"
  if [ -e $nf ] ; then
    
    if [ `echo $bn | fgrep "pwd"` ] ; then
        testpwd=1
        moebatch -exec "run ['pwdbatch.svl', ['$nf', 'QMScore', 15]]" -exit
    else
        testpwd=0
    fi

    moebatch -exec "run ['qbmdbdiff.svl', ['qbmdbdiff-error.log', '$bf', '$nf', 'all', 0.5, $testpwd]]" -exit

  else
    echo "$nf doesn't exist."
  fi

done



