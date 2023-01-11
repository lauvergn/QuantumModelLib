#!/bin/bash

EXTLIB_TYPE=$1

echo "In get_dnSVM.sh"
pwd


SAVE_version=Save_AD_dnSVM-3.2
LOC_version=AD_dnSVM


rm -rf AD_dnSVM* #always remove the link


#latest release
#latest HEAD version (dev version)
 version=https://github.com/lauvergn/AD_dnSVM/archive/refs/tags/v3.2.zip


test -z $EXTLIB_TYPE       &&    curl -LJ $version --output $LOC_version.zip
test $EXTLIB_TYPE != 'loc' &&    curl -LJ $version --output $LOC_version.zip

test -e $LOC_version.zip && echo $LOC_version.zip file exist || cp $SAVE_version.zip $LOC_version.zip

unzip $LOC_version.zip
rm -f $LOC_version.zip


LIBDIR=`ls -d AD_dnSVM*`
#echo $LIBDIR

ln -s $LIBDIR $LOC_version

echo "End get_dnSVM.sh"

