#!/usr/bin/env bash

case $TRAVIS_OS_NAME in
    osx)
        brew update
        for pkg in metis bash gcc; do
            if [ x$(brew list | fgrep $pkg) != x$pkg ]; then
                brew install $pkg
            fi
        done
        ;;
    windows)
        [[ ! -f C:/tools/msys64/msys2_shell.cmd ]] && rm -rf C:/tools/msys64
        choco uninstall -y mingw
        choco upgrade --no-progress -y msys2
        $BASH pacman --sync --noconfirm --needed mingw-w64-x86_64-toolchain
        $BASH pacman -S mingw-w64-x86_64-lapack --noconfirm
        $BASH pacman -S mingw-w64-x86_64-winpthreads-git --noconfirm
        $BASH pacman -S mingw-w64-x86_64-readline --noconfirm
        $BASH pacman -S mingw-w64-x86_64-suitesparse --noconfirm
        $BASH pacman -S mingw-w64-x86_64-metis --noconfirm
        $BASH pacman -S make wget tar patch dos2unix diffutils --noconfirm
        $BASH pacman -S pkg-config git zip unzip --noconfirm
        taskkill //IM gpg-agent.exe //F  # https://travis-ci.community/t/4967
        ;;
esac

