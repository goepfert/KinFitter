#!/bin/sh

# this needs to be sourced
# set -x

test $# -eq 1 && test $1 == "-h" && echo "Usage: setANALpath [source dir] [BFARCH]" && return 1;

if [ -z "$ROOTSYS" ]
then
    echo "No ROOT version is set. Use 'setROOTpath' before running setANALpath"
    return 1;
fi

export ROOTRELEASE=`root-config --version | sed -e 's/\//\./g'`

# check the src dir $1

if [  $# -gt 0 ] && [ ! -d $1 ]
then
    echo "Directory $1 doesn't exist!"
    return 3;
fi

if [ -z "$1" ]; then
    export SRC=`pwd`
else
    export SRC=$1
fi

if [ $# -eq 2 ]; then
    export BFARCH=$2
else
    if [ -z "$BFARCH" ]; then
	#try to guess the BFARCH
	kernel=`uname -r|awk -F"." '{print $1$2}'`

	which lsb_release 1> /dev/null 2> /dev/null; excode=$?

	if  [ $excode -ne 0 ]; then
	    echo "Cannot guess BFARCH because lsb_release cannot be found! Please specify a BFARCH on the command line!"
	    return 2;
	fi

	if lsb_release -d | grep -q Scientific ; then
	    vendor=SL
	    release=`lsb_release -r |awk '{print $2}' |awk -F"." '{print $1}'`
	fi
	if lsb_release -d | grep -q RedHat ; then
	    vendor=RHEL
	    release=`lsb_release -r |awk '{print $2}'`
	fi

	if uname -m |grep -q -E i?86 ; then
	    arch=i386
	fi
	if uname -m |grep -q 64 ; then
	    arch=x86_64
	fi

	gccver=`gcc -dumpversion| sed -e 's/\.//g'`

	placeholder1="_"
	placeholder2="_gcc"
	BFARCH="Linux$kernel$vendor$release$placeholder1$arch$placeholder2$gccver"
	export BFARCH
    fi
fi

if [ -n $LIB ]; then
    oldlib=$LIB
fi


export LIB=`echo $SRC/lib/$BFARCH/$ROOTRELEASE |sed -e 's%\/\/%\/%g'`
if [ x$oldlib != x ]; then
    sedlib=`echo $oldlib | sed  -e 's/\_/\\\_/g'  -e 's%\/%\\\/%g' -e 's%\.%\\\.%g'`
    sedstring='s%:[^:]*'$sedlib'[^:]*:%:%g'
    export LD_LIBRARY_PATH=`echo :$LD_LIBRARY_PATH: | sed -e $sedstring  -e 's/^://g' -e 's% .$%%g' -e 's%:$%%g' -e 's%\/\/%\/%g'`:$LIB
fi

if [ -z $oldlib ]; then
    export LD_LIBRARY_PATH=`echo :$LD_LIBRARY_PATH: | sed  -e 's/^://g' -e 's% .$%%g' -e 's%:$%%g' -e 's%\/\/%\/%g'`:$LIB
fi


echo "Setup build environment for ROOT version $ROOTRELEASE in $SRC with BFARCH $BFARCH !"
return 0
