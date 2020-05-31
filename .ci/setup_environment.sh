#!/usr/bin/env bash

case $TRAVIS_OS_NAME in
    osx)
        brew update
        brew install metis bash gcc
        if [ $CC = gcc ]; then
            export CC=gcc-9
            export CXX=g++-9
        fi
        ;;
    windows)
        [[ ! -f C:/tools/msys64/msys2_shell.cmd ]] && rm -rf C:/tools/msys64
        choco uninstall -y mingw
        choco upgrade --no-progress -y msys2
        export msys2='cmd //C RefreshEnv.cmd '
        export msys2+='& set MSYS=winsymlinks:nativestrict '
        export msys2+='& C:\\tools\\msys64\\msys2_shell.cmd -defterm -no-start'
        export mingw64="$msys2 -mingw64 -full-path -here -c "\"\$@"\" --"
        export msys2+=" -msys2 -c "\"\$@"\" --"
        $msys2 pacman --sync --noconfirm --needed mingw-w64-x86_64-toolchain
        $msys2 pacman -S mingw-w64-x86_64-lapack --noconfirm
        $msys2 pacman -S mingw-w64-x86_64-winpthreads-git --noconfirm
        $msys2 pacman -S mingw-w64-x86_64-readline --noconfirm
        $msys2 pacman -S mingw-w64-x86_64-suitesparse --noconfirm
        $msys2 pacman -S mingw-w64-x86_64-metis --noconfirm
        $msys2 pacman -S make wget tar patch dos2unix diffutils --noconfirm
        $msys2 pacman -S pkg-config git zip unzip --noconfirm
        taskkill //IM gpg-agent.exe //F  # https://travis-ci.community/t/4967
        export PATH=/C/tools/msys64/mingw64/bin:$PATH
        export MAKE=mingw32-make  # so that Autotools can find it
        ;;
esac
case $CC in
    gcc*)
        export CCVERSION=$($CC -dumpversion)
        ;;
    clang)
        export CCVERSION=$(clang --version |& fgrep version |& \
                           sed "s/.*version \([0-9]*\.[0-9]*\).*/\1/")
        ;;
esac
declare -a DBG_ARGS
export DBG_ARGS=()
if [ $DEBUG = "true" ]; then
    export DBG_ARGS=( --enable-debug )
    export DBGN="-dbg"
    export CXXFLAGS=( -Og -g)
fi
if [ $ASAN = "true" ]; then
    export ASN="-asan"
    export CXXFLAGS="${CXXFLAGS} -fsanitize=address"
    export LDFLAGS="-lasan"
fi
declare -a DBG_ARGS
export ADD_ARGS=()
if [ $BUILD_STATIC = "true" ]; then
    export STATIC=( -static )
    export ADD_ARGS=( --static --with-lapack='-llapack -lblas -lgfortran -lquadmath -lm' )
fi        
declare -a COMMON_ARGS
export COMMON_ARGS=( --no-prompt --verbosity 2 --tests main --enable-relocatable )
export PLATFORM=$TRAVIS_OS_NAME${OSX:-}-x86_64-$CC$CCVERSION
export PROJECT_URL=https://github.com/$TRAVIS_REPO_SLUG
