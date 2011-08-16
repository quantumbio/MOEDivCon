#!/bin/bash







# where the HDFView is installed
if [ ! "$QBHOME" ]; then
    PRG=$0
    
    progdir=`dirname "$PRG"`
    if [ "$progdir" = "." ]; then
	    progdir=`pwd`
    fi
    
    if [ "`echo "$progdir" | grep -s '^/'`" ]; then
        progdir=$progdir
    else
        progdir=`pwd`/"$progdir"
    fi
    
    QBHOME=`echo "$progdir" | sed 's:/bin$::'`
    export QBHOME
fi

if [ -z "$QBENVSET" ]; then
    . $QBHOME/etc/qbenv.sh
fi

HDFVIEW_HOME=$QBHOME
export HDFVIEW_HOME

# where Java is installed (requires jdk1.4.x or above), e.g. /usr/jdk1.4.2/bin
java_bin=`which java`
if [ ! -n "$java_bin" ]; then
    echo "ERROR! Java $JAVAVERSIONREQ not installed."
    exit 1
fi
java_version=`java -version 2>&1 | head -1 | awk -F "\"" '{print $2}' | awk -F "." '{print $1"."$2}'`
if [ $java_version -lt $JAVAVERSIONREQ ]; then
    echo "ERROR! Java $JAVAVERSIONREQ required."
    exit 1
fi
JAVAPATH=`dirname "$java_bin"`
export JAVAPATH

CPATH=$HDFVIEW_HOME/$PLTFM_VER/HDFView.jar
CPATH=$CPATH:$HDFVIEW_HOME/$PLTFM_VER/lib/netcdf.jar:$HDFVIEW_HOME/$PLTFM_VER/lib/fits.jar


TEST=/usr/bin/test
if [ ! -x /usr/bin/test ] 
then
TEST=`which test`
fi

if $TEST -z "$CLASSPATH"; then
	CLASSPATH=""
fi
CLASSPATH=$CPATH":"$CLASSPATH
export CLASSPATH

if $TEST -n "$JAVAPATH" ; then
	PATH=$JAVAPATH":"$PATH
	export PATH
fi

$JAVAPATH/java -Xmx1000m -Djava.library.path=$LD_LIBRARY_PATH ncsa.hdf.view.HDFView -root $HDFVIEW_HOME $*
