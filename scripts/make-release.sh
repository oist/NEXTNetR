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

echo "Extracting archive to .buildpkg/" >&2
rm -rf .buildpkg
mkdir .buildpkg
(cd .buildpkg; tar xzf ../archives/NEXTNetR-v$ver-full.tar.gz)

echo "Building proper R package archives/NEXTNetR-v$ver-pkg.tar.gz" >&2
R CMD build .buildpkg/NEXTNetR
mv NEXTNetR_$ver.tar.gz archives/NEXTNetR-v$ver-pkg.tar.gz
rm -r .buildpkg
