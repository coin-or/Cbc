#!/usr/bin/env bash

case $TRAVIS_OS_NAME in
    osx)
        brew update
        brew install metis bash gcc
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
        ;;
esac

