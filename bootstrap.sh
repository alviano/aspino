#!/bin/sh

DIRNAME=`dirname $0`
OLDPATH=`pwd`

echo "Removing directory src/glucose-syrup (if it exists)"
rm -rf $DIRNAME/src/glucose-syrup

echo "Downloading glucose 4.0..."
if [ ! -e $DIRNAME/patches/glucose-syrup.tgz ]; then
    wget -O $DIRNAME/patches/glucose-syrup.tgz http://www.labri.fr/perso/lsimon/downloads/softwares/glucose-syrup.tgz
fi

echo "Patching glucose 4.0..."
cd $DIRNAME/patches; tar xf glucose-syrup.tgz; cd $OLDPATH
cd $DIRNAME/patches/glucose-syrup; patch -p1 <../glucose-syrup.4.0.patch.1; cd $OLDPATH
mv $DIRNAME/patches/glucose-syrup $DIRNAME/src/

echo "Downloading gflags 2.1.1..."
mkdir -p lib
if [ ! -e $DIRNAME/lib/gflags-2.1.1.tar.gz ]; then
    wget -O $DIRNAME/lib/gflags-2.1.1.tar.gz https://github.com/schuhschuh/gflags/archive/v2.1.1.tar.gz
fi

echo "Building gflags 2.1.1..."
cd $DIRNAME/lib; tar xf gflags-2.1.1.tar.gz; cd $OLDPATH
cd $DIRNAME/lib/gflags-2.1.1; mkdir -p build; cd build; cmake ..; make; cd $OLDPATH

echo
echo "You can now run make!"
