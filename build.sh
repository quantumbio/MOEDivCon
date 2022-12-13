#!/bin/bash

function pause(){
   read -n 1 -p "Press any key to continue..."
}

export JAVA_HOME=/home/roger/jdk-11.0.16.1
export ANT_HOME=/home/roger/apache/apache-ant-1.10.12
export METRO_HOME=/home/roger/java.net/metro
export PATH=$JAVA_HOME/bin:$ANT_HOME/bin:$METRO_HOME/bin:/opt/gcc/current/bin:/home/roger/Software/cmake-2.8.10.2/dist/bin:$PATH
#export LD_LIBRARY_PATH=/share/apps/intel/linux/composer_xe_2011_sp1.9.293/mkl/lib/intel64:/share/apps/intel/linux/composer_xe_2011_sp1.9.293/compiler/lib/intel64:/home/roger/NetBeansProjects/openmpi/lib:/home/roger/NetBeansProjects/antlr/lib:/home/roger/NetBeansProjects/boost/lib:/home/roger/NetBeansProjects/lmx/lib:/home/roger/NetBeansProjects/hdf5/lib:/home/roger/NetBeansProjects/log4cplus/lib:/home/roger/NetBeansProjects/icu4c/lib:/home/roger/NetBeansProjects/zlib/lib:/home/roger/NetBeansProjects/bzip2/lib:/home/roger/NetBeansProjects/FockianIntegrals/dist/Release/GNU-Linux-x86:/home/roger/NetBeansProjects/OOBackbone/dist/Release/GNU-Linux-x86:/home/roger/NetBeansProjects/JPRoothaan/dist/Release/GNU-Linux-x86:/home/roger/NetBeansProjects/QBPersistence/dist/Release/GNU-Linux-x86:/opt/gcc/current/lib64:/lib64:$LD_LIBRARY_PATH

export libs_jaxb_classpath="`pwd`/lib/jaxb-impl.jar:`pwd`/lib/jaxb-xjc.jar:`pwd`/lib/jaxb1-impl.jar:`pwd`/lib/activation.jar:`pwd`/lib/jaxb-api.jar"
/home/roger/apache/apache-ant-1.10.12/bin/ant clean
/home/roger/apache/apache-ant-1.10.12/bin/ant jar

pause
