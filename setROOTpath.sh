#!/bin/sh
# find PAF LD_LIBRARY_PATH and replace by specified PATH.
# usage: setROOTpath /path/to/root

if [ -d "$1" ]; then
    export PATH=`echo :$PATH:| sed -e 's/:[^:]*root[^:]*:/:/g' -e 's/^://g' -e 's/ .$//g' -e 's/:$//g'`:$1/bin
    export LD_LIBRARY_PATH=`echo :$LD_LIBRARY_PATH:| sed -e 's/:[^:]*root[^:]*:/:/g' -e 's/^://g' -e 's/ .$//g' -e 's/:$//g'`:$1/lib
    export ROOTSYS=$1
    export BBRROOTSYS=$ROOTSYS
else
    echo Directory does not exist: $1
fi
