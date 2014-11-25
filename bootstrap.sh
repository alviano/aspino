#!/bin/sh

DIRNAME=`dirname $0`
OLDPATH=`pwd`

echo "Removing directory src/glucose-syrup (if it exists)"
rm -rf $DIRNAME/src/glucose-syrup

echo "Downloading glucose 4.0..."
wget -O $DIRNAME/patches/glucose-syrup.tgz http://www.labri.fr/perso/lsimon/downloads/softwares/glucose-syrup.tgz

echo "Patching glucose 4.0..."
cd $DIRNAME/patches; tar xf glucose-syrup.tgz; cd $OLDPATH
cd $DIRNAME/patches/glucose-syrup; patch -p1 <../glucose-syrup.4.0.patch; cd $OLDPATH
mv $DIRNAME/patches/glucose-syrup $DIRNAME/src/
rm -f $DIRNAME/patches/glucose-syrup.tgz

echo
echo "You can now run make!"
