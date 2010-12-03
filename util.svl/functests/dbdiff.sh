#!/bin/bash
#

if [ $# != 2 ] ; then
  echo "Usage:  dbdiff.sh basedir newdir"
  exit
fi

bfiles=`/bin/ls -1 $1/*.mdb`

for bf in $bfiles ; do
  bn=`basename $bf`
  nf="${2}/$bn"
  if [ -e $nf ] ; then
    
    if [ `echo $bn | fgrep "pwd"` ] ; then
        testpwd=1
    else
        testpwd=0
    fi

    moebatch -exec "run ['qbmdbdiff.svl', ['qbmdbdiff-error.log', '$bf', '$nf', 'all', 0.25, $testpwd]]" -exit

  else
    echo "$nf doesn't exist."
  fi

done



