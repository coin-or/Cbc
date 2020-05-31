#!/usr/bin/env bash

if [ $TRAVIS_BRANCH = "master" ]; then 
    VERSION=master 
else 
    VERSION=`echo $TRAVIS_BRANCH | cut -d "/" -f 2`
fi

case $TRAVIS_OS_NAME in
    linux)
        TGZ_FILE=$PROJECT-$VERSION${STATIC:-}-$TRAVIS_DIST-$PLATFORM${DBGN:-}${ASN:-}.tgz
        ;;
    *)
        TGZ_FILE=$PROJECT-$VERSION${STATIC:-}-$PLATFORM${DBGN:-}${ASN:-}.tgz
esac
echo $TGZ_FILE
cp $PROJECT/README.md $PROJECT/LICENSE $PROJECT/AUTHORS dist
cd dist
if [ $TRAVIS_OS_NAME = "linux" ] && [ $TRAVIS_DIST = precise ]; then
    # adding required libraries
    cp /usr/lib/lapack/liblapack.so.3gf lib/
    cp /usr/lib/x86_64-linux-gnu/libgfortran.so.5 lib/
    cp /usr/lib/x86_64-linux-gnu/libgfortran.so.3 lib/
    cp /usr/lib/x86_64-linux-gnu/libstdc++.so.6 lib/
    cp /usr/lib/libblas/libblas.so.3gf lib/
    cp /lib/x86_64-linux-gnu/libreadline.so.6 lib/
    cp /lib/x86_64-linux-gnu/libbz2.so.1.0 lib/
    cp /lib/x86_64-linux-gnu/libtinfo.so.5 lib/
    cp /usr/lib/x86_64-linux-gnu/libquadmath.so.0 lib/
    cp /usr/lib/libcholmod.so.1.7.1 lib/
    cp /usr/lib/libamd.so.2.2.0 lib/
    cp /usr/lib/libcolamd.so.2.7.1 lib/
    cp /usr/lib/libmetis-edf.so.* lib/
    
    for libfile in lib/*.so*
    do
        chrpath -r ./ $libfile
    done
    
    chrpath -r \$\ORIGIN/../lib/ bin/cbc
    chrpath -r \$\ORIGIN/../lib/ bin/clp
    chrpath -r \$\ORIGIN/../lib/ bin/glpsol
fi
rm lib/*.la
echo $TGZ_FILE
set -x
tar -czvf $TGZ_FILE lib/* bin/* include/* share/* README.md LICENSE AUTHORS
curl -T $TGZ_FILE -utkralphs:$BINTRAY_API -H "X-Bintray-Publish:1" \
     -H "X-Bintray-Override:1" \
     https://api.bintray.com/content/coin-or/download/$PROJECT/$VERSION/$TGZ_FILE
while [[ ! $(curl -utkralphs:$BINTRAY_API -X PUT \
                 -H "Content-Type: application/json" \
                 -d'{"list_in_downloads":true}' \
                 https://api.bintray.com/file_metadata/coin-or/download/$TGZ_FILE) =~ success ]];
do 
    echo "Command failed...trying again 10 seconds"
    sleep 10
done
cd ..
