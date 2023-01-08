#!/bin/bash

echo "In get_dnSVM.sh"
pwd


FC=$1
OPT=$2
OMP=$3
LAPACK=$4
ExtLibDIR=$5

ext_obj="_"$FC"_opt"$OPT"_omp"$OMP

SAVE_version=Save_AD_dnSVM-AD_dnSVM_dev
LOC_version=AD_dnSVM

test -f $ExtLibDIR/$LOC_version/"libAD_dnSVM"$ext_obj.a && exit 0


rm -rf AD_dnSVM* #always remove the link


#latest release
#latest HEAD version (dev version)
 version=https://github.com/lauvergn/AD_dnSVM/archive/refs/heads/AD_dnSVM_dev.zip


curl -LJ $version --output $LOC_version.zip
test -e $LOC_version.zip && echo $LOC_version.zip file exist || cp $SAVE_version.zip $LOC_version.zip

unzip $LOC_version.zip
rm -f $LOC_version.zip

LIBDIR=`ls -d AD_dnSVM*`
#echo $LIBDIR

ln -s $LIBDIR $LOC_version

cd $LIBDIR
  make lib FC=$FC OPT=$OPT OMP=$OMP LAPACK=$LAPACK ExtLibDIR=$ExtLibDIR

echo "End get_dnSVM.sh"

