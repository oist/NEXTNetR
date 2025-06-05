#!/bin/bash
set -e
set -o pipefail 

ver=$1
if [ "$ver" == "" ]; then
	echo "Usage: $0 version" >&2
	exit 1
fi

if [ $(git ls-remote --tags origin v$ver | wc -l) != 0 ]; then
	echo "Re-using exiting tag for version $ver" >&2
else
	./scripts/tag-release.sh $ver
fi

echo "Building full archive for version $ver (including submodules)" >&2
mkdir -p archives
rm -f  archives/NEXTNetR-v$ver-full.tar
rm -f  archives/NEXTNetR-v$ver-full.tar.gz
./scripts/git-archive-all.sh --tree-ish v$ver --prefix NEXTNetR/ archives/NEXTNetR-v$ver-full.tar
gzip --best archives/NEXTNetR-v$ver-full.tar

echo "Building source archive for version $ver" >&2
rm -rf .buildpkg
mkdir .buildpkg
cd .buildpkg
tar xzf ../archives/NEXTNetR-v$ver-full.tar.gz
tar czf ../archives/NEXTNetR-v$ver-src.tar.gz \
	NEXTNetR/DESCRIPTION \
	NEXTNetR/NAMESPACE \
	NEXTNetR/configure \
	NEXTNetR/src/Makevars* \
	NEXTNetR/src/*.{cpp,h,hpp} \
	NEXTNetR/ext/nextnet/nextnet/*.{cpp,h} \
	NEXTNetR/ext/nextnet/nextnet/pstream/*.{cpp,h} \
	NEXTNetR/ext/nextnet/ext/prio_queue/prio_queue.hpp \
	NEXTNetR/ext/nextnet/ext/dyndist/dyndist/*.h \
	NEXTNetR/R/*.R \
	NEXTNetR/man/*.Rd \
	NEXTNetR/examples/*.R \
	NEXTNetR/inst/CITATION \
	NEXTNetR/tests/** \
	NEXTNetR/vignettes/*.{Rmd,nw} \
	NEXTNetR/.Rbuildignore \
	NEXTNetR/LICENSE.md \
	NEXTNetR/README.md
cd ../
rm -rf .buildpkg


echo "Building proper R package archives/NEXTNetR-v$ver-pkg.tar.gz" >&2
rm -rf .buildpkg
mkdir .buildpkg
cd .buildpkg
tar xzf ../archives/NEXTNetR-v$ver-src.tar.gz
R CMD BUILD NEXTNetR
mv NEXTNetR_$ver.tar.gz ../archives/NEXTNetR-v$ver-pkg.tar.gz
cd ../
rm -r .buildpkg
