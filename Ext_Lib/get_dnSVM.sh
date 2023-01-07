rm -rf AD_dnSVM* #always remove the link

SAVE_Version=Save_AD_dnSVM-AD_dnSVM_dev

LOC_version=AD_dnSVM-loc

#latest release
#latest HEAD version
 version=https://github.com/lauvergn/AD_dnSVM/archive/refs/heads/main.zip


#curl -LJ $version --output $LOC_version.zip
test -e $LOC_version.zip && echo $LOC_version.zip file exist || cp $SAVE_Version.zip $LOC_version.zip

unzip $LOC_version.zip
rm -f $LOC_version.zip

LIBDIR=`ls -d AD_dnSVM*`
#echo $LIBDIR

ln -s $LIBDIR $LOC_version

cd $LIBDIR
  pwd
  make lib OPT=$1 OMP=$2 LAPACK=$3
