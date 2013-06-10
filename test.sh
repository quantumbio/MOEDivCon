#!/bin/bash
function pause(){
   read -n 1 -p "Press any key to continue..."
}

date

export PROJECT_PATH=/home/roger/DivConDiscoverySuite-b1251
export JAVA_HOME=/home/roger/jdk1.7.0_13
export ANT_HOME=/home/roger/apache/apache-ant-1.8.4
export OOBACKBONE_HOME=/home/roger/NetBeansProjects/OOBackbone
export PERSISTENCE_HOME=/home/roger/NetBeansProjects/QBPersistence
export QMECHANIC_HOME=/home/roger/testing-grounds/qmechanic

export PATH=/share/apps/MOE/MOE-2012/bin:$JAVA_HOME/bin:$ANT_HOME/bin:/home/roger/NetBeansProjects/openmpi/bin:$QMECHANIC_HOME/bin:$PROJECT_PATH/bin:$PATH
export LD_LIBRARY_PATH=/home/roger/Software/hdf-java-2.9/lib/linux:$QMECHANIC_HOME/linux-x86_64/lib:/share/apps/intel/linux/composer_xe_2013.2.146/mkl/lib/intel64:/share/apps/intel/linux/composer_xe_2013.2.146/compiler/lib/intel64:/home/roger/NetBeansProjects/openmpi/lib:/home/roger/NetBeansProjects/antlr/lib:/home/roger/NetBeansProjects/boost/lib:/home/roger/NetBeansProjects/lmx/lib:/home/roger/NetBeansProjects/hdf5/lib:/home/roger/NetBeansProjects/log4cplus/lib:/home/roger/NetBeansProjects/icu4c/lib:/home/roger/NetBeansProjects/zlib/lib:/home/roger/NetBeansProjects/bzip2/lib:/home/roger/NetBeansProjects/FockianIntegrals/dist/Release/GNU-Linux-x86:/home/roger/NetBeansProjects/OOBackbone/dist/Release/GNU-Linux-x86:/home/roger/NetBeansProjects/JPRoothaan/dist/Release/GNU-Linux-x86:/home/roger/NetBeansProjects/QBPersistence/dist/Release/GNU-Linux-x86:/home/roger/NetBeansProjects/JNIDivcon/native/JNIDivcon/dist/Release/GNU-Linux-x86:/opt/gcc/current/lib64:/lib64:$LD_LIBRARY_PATH

export CLASSPATH=`pwd`/dist/MOEDivcon.jar:$QMECHANIC_HOME/linux-x86_64/HDFView.jar:$QMECHANIC_HOME/linux-x86_64/QBReporter.jar:/share/apps/MOE/MOE-2012/java/svljava.jar:$CLASSPATH
export SVL_CLASSPATH=$CLASSPATH
export OPAL_PREFIX=$QMECHANIC_HOME/linux-x86_64
#ant test-single -DOPAL_PREFIX=$OPAL_PREFIX -DQMECHANIC_HOME=$QMECHANIC_HOME -Djava.library.path=$LD_LIBRARY_PATH -Dtest.includes=HDF5Tests.java -Djavac.includes=HDF5Tests.java -Dtest.class=HDF5Tests

rm roothaan.h5
export MOE_SVL_LOAD=/home/roger/NetBeansProjects/MOEDivcon/svl
#qbmoebatch -exec "run ['svl/hdf5interface.svl', ['/home/roger/testing-grounds/s59/Moe interface/s59.pm6.sp.xml.h5']]" -exit
#qbmoebatch -exec "qbTesth5Main['/home/roger/testing-grounds/s59/Moe interface/s59.pm6.sp.xml.h5']" -exit
cp /home/roger/NetBeansProjects/MOEDivcon/dist/MOEDivcon.jar /home/roger/DivConDiscoverySuite-b1251/linux-x86_64
qbmoebatch -exec "qbTestPDBtoH5['./roothaan.h5', '/home/roger/NetBeansProjects/OOBackbone/tests/examples/3FVA.pdb']" -exit

#qbmoe

pause