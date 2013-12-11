#!/bin/bash
function pause(){
   read -n 1 -p "Press any key to continue..."
}

date

export PROJECT_PATH=/home/roger/DivConDiscoverySuite-b1761
#export JAVA_HOME=/home/roger/jdk1.7.0_25
#export ANT_HOME=/home/roger/apache/apache-ant-1.8.4
#export OOBACKBONE_HOME=/home/roger/NetBeansProjects/OOBackbone
#export PERSISTENCE_HOME=/home/roger/NetBeansProjects/QBPersistence
#export QMECHANIC_HOME=/home/roger/DivConDiscoverySuite-b1418

export PATH=$PROJECT_PATH/bin:$PATH
#qbmoebatch -exec "qbTesth5Main['/home/roger/testing-grounds/Lance/1LRI/1LRI-NMR.h5']" -exit

export LD_LIBRARY_PATH=/home/roger/Software/hdf-java-2.9/lib/linux:$QMECHANIC_HOME/linux-x86_64/lib:/share/apps/intel/linux/composer_xe_2013.5.192/mkl/lib/intel64:/share/apps/intel/linux/composer_xe_2013.5.192/compiler/lib/intel64:/home/roger/NetBeansProjects/openmpi/lib:/home/roger/NetBeansProjects/antlr/lib:/home/roger/NetBeansProjects/boost/lib:/home/roger/NetBeansProjects/lmx/lib:/home/roger/NetBeansProjects/hdf5/lib:/home/roger/NetBeansProjects/log4cplus/lib:/home/roger/NetBeansProjects/icu4c/lib:/home/roger/NetBeansProjects/zlib/lib:/home/roger/NetBeansProjects/bzip2/lib:/home/roger/NetBeansProjects/FockianIntegrals/dist/Release/GNU-Linux-x86:/home/roger/NetBeansProjects/OOBackbone/dist/Release/GNU-Linux-x86:/home/roger/NetBeansProjects/JPRoothaan/dist/Release/GNU-Linux-x86:/home/roger/NetBeansProjects/QBPersistence/dist/Release/GNU-Linux-x86:/home/roger/NetBeansProjects/JNIDivcon/native/JNIDivcon/dist/Release/GNU-Linux-x86:/opt/gcc/current/lib64:/lib64:$LD_LIBRARY_PATH
#export LD_LIBRARY_PATH=/home/roger/NetBeansProjects/icu4c/lib:/opt/gcc/current/lib64:/lib64:$LD_LIBRARY_PATH
#export CLASSPATH=`pwd`/dist/MOEDivcon.jar:$QMECHANIC_HOME/linux-x86_64/HDFView.jar:$QMECHANIC_HOME/linux-x86_64/QBReporter.jar:/share/apps/MOE/MOE-2012/java/svljava.jar:$CLASSPATH
#export SVL_CLASSPATH=$CLASSPATH
#export OPAL_PREFIX=$QMECHANIC_HOME/linux-x86_64
#ant test-single -DOPAL_PREFIX=$OPAL_PREFIX -DQMECHANIC_HOME=$QMECHANIC_HOME -Djava.library.path=$LD_LIBRARY_PATH -Dtest.includes=HDF5Tests.java -Djavac.includes=HDF5Tests.java -Dtest.class=HDF5Tests

rm roothaan.h5
export MOE_SVL_LOAD=/home/roger/NetBeansProjects/MOEDivcon/svl
cp /home/roger/NetBeansProjects/MOEDivcon/dist/MOEDivcon.jar $PROJECT_PATH/linux-x86_64
cp /home/roger/NetBeansProjects/MOEDivcon/svl/hdf5interface.svl $PROJECT_PATH/svl
cp /home/roger/NetBeansProjects/JNIDivcon/dist/JNIDivcon.jar $PROJECT_PATH/linux-x86_64/lib

#qbmoebatch -exec "qbTesth5Main['/home/roger/testing-grounds/Josh/abbott/x_mcl1_a1107644/x_mcl1_a1107644_protonated.pdb.h5']" -exit
#qbmoebatch -exec "qbTesth5Main['/home/roger/testing-grounds/s59/Moe interface/s59.pm6.sp.xml.h5']" -exit
#qbmoebatch -exec "qbTesth5Main['/home/roger/testing-grounds/Josh/s59/s59.pm6.sp.xml.h5']" -exit
#qbmoebatch -exec "qbTestPDBtoH5['./roothaan.h5', '/home/roger/NetBeansProjects/OOBackbone/tests/examples/3FVA.pdb']" -exit
#qbmoebatch -exec "qbTesth5Main['/home/roger/testing-grounds/Cl4Y2/roothaan.h5']" -exit
#qbmoebatch -exec "qbTesth5Main['/home/roger/testing-grounds/Josh/CSPtest/shadowfax/2YXJ_out.xml.h5']" -exit
#qbmoebatch -exec "qbTesth5Main['/home/roger/testing-grounds/Josh/CSPtest/frodo/roothaan.h5']" -exit
#qbmoebatch -exec "qbTesth5Main['/home/roger/testing-grounds/Lance/2YXJ/Tut#1/targ_CSP/2YXJ_target_shifts_dc.h5']" -exit
#qbmoebatch -exec "qbTesth5Main['/home/lance/tmp/targ_NMRScore.node/2YXJ_target_shifts_dc.h5']" -exit
qbmoebatch -exec "qbTesth5Main['./2YXJ_target_shifts_dc.h5']" -exit

#qbmoe

pause
