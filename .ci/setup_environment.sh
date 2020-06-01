#!/usr/bin/env bash

case $CC in
    gcc*)
        if [ $TRAVIS_OS_NAME = osx ]; then
            export CC=gcc-9
            export CXX=g++-9
            CCVERSION=9
        else
            export CCVERSION=$($CC -dumpversion)
        fi
        ;;
    clang)
        export CCVERSION=$(clang --version | fgrep version | \
                           sed "s/.*version \([0-9]*\.[0-9]*\).*/\1/")
        ;;
esac
declare -a DBG_ARGS
export DBG_ARGS=()
if [ x$DEBUG = x"true" ]; then
    export DBG_ARGS=( --enable-debug )
    export DBGN="-dbg"
    export CXXFLAGS=( -Og -g)
fi
if [ x$ASAN = x"true" ]; then
    export ASN="-asan"
    export CXXFLAGS="${CXXFLAGS} -fsanitize=address"
    export LDFLAGS="-lasan"
fi
declare -a ADD_ARGS
export ADD_ARGS=()
if [ x$BUILD_STATIC = x"true" ]; then
    export STATIC=( -static )
    export ADD_ARGS=( --static --with-lapack='-llapack -lblas -lgfortran -lquadmath -lm' )
fi        
declare -a COMMON_ARGS
export COMMON_ARGS=( --no-prompt --verbosity 2 --tests main --enable-relocatable )
export PLATFORM=$TRAVIS_OS_NAME${OSX:-}-x86_64-$CC$CCVERSION
export PROJECT_URL=https://github.com/$TRAVIS_REPO_SLUG
if [ $TRAVIS_OS_NAME = windows ]; then
    export PATH=/C/tools/msys64/mingw64/bin:$PATH
fi
