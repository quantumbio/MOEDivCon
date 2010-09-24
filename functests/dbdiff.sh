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
    moebatch -exec "run ['qbmdbdiff.svl', ['qbmdbdiff-error.log', '$bf', '$nf', 'all', 0.05]]" -exit
  else
    echo "$nf doesn't exist."
  fi

done



